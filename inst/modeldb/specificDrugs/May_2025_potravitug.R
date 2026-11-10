May_2025_potravitug <- function() {
  description <- paste(
    "QSP. Semi-mechanistic PK/PD model of potravitug, a fully human IgG1",
    "against the BK polyomavirus (BKPyV) major capsid protein VP1, in",
    "healthy adults (May 2025). Three drug compartments (central blood,",
    "peripheral tissues, kidney interstitial space) with first-order",
    "catabolism in every compartment and Pdist/Tdist inter-compartmental",
    "exchange; VP1 (BKPyV titer x 360 VP1 per virion) turns over in blood",
    "and kidney interstitial space and binds potravitug 1:1, the complex",
    "being cleared at the VP1 rate. Deterministic typical-subject (70 kg)",
    "model calibrated on cohort-mean SAD data; no IIV or residual error.",
    sep = " "
  )
  reference <- paste(
    "May M, Bialek-Waldmann JK, Wright A, Gruver S, Grant J, Eicher B,",
    "Seidenberg J, Parzmair GP, Beck J. First-In-Human, Randomized,",
    "Placebo-Controlled, SAD and MAD Trial to Evaluate Safety, Tolerability,",
    "and PK/PD Modeling of Potravitug in Healthy Adults. Clin Transl Sci.",
    "2025;18(8):e70316. doi:10.1111/cts.70316. Model parameters from",
    "Appendix S1 Table S1.",
    sep = " "
  )
  vignette <- "May_2025_potravitug"
  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "ug/mL"
  )

  # VP1 free and potravitug-bound states in the kidney interstitial space.
  # The bare `target` / `complex` carry the same species in blood.
  paper_specific_compartments <- c("target_kidney", "complex_kidney")

  compartmentData <- list(
    central = list(analyte = "potravitug (free)", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "potravitug (free)", units = "mg", specimen = "tissue", verified = TRUE),
    kidney = list(
      analyte = "potravitug (free), kidney interstitial space",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    ),
    target = list(analyte = "BKPyV VP1 (free)", units = "nmol", specimen = "whole blood", verified = TRUE),
    complex = list(analyte = "potravitug-VP1 complex", units = "nmol", specimen = "whole blood", verified = TRUE),
    target_kidney = list(
      analyte = "BKPyV VP1 (free), kidney interstitial space",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    complex_kidney = list(
      analyte = "potravitug-VP1 complex, kidney interstitial space",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    )
  )

  covariateData <- list()

  population <- list(
    species = "human",
    n_subjects = 40L,
    n_studies = 1L,
    age_range = "18-50 years",
    weight_range = "51.2-100.5 kg",
    sex_female_pct = 50,
    race_ethnicity = c(White = 82.5, Asian = 10, Black = 2.5, `Pacific Islander` = 2.5, `Asian and White` = 2.5),
    disease_state = "Healthy adult volunteers (BMI 18-32 kg/m^2); the VP1 sub-model represents an average BKPyV titer as in kidney-transplant recipients.",
    dose_range = "Single IV infusion 100, 500, 1000 or 2000 mg (SAD, n = 12 active); 500, 1000 or 2000 mg IV every 4 weeks x 4 (MAD, n = 18 active).",
    regions = "Australia (single centre, Adelaide)",
    notes = paste(
      "First-in-human trial NCT05358106; demographics from May 2025 Table 1",
      "(40 randomized, of whom 30 received potravitug). The PK/PD model was",
      "calibrated on the SAD cohort-mean serum concentrations for a typical",
      "70 kg subject and validated against the MAD cohort means; it has no",
      "covariates and no between-subject variability."
    )
  )

  ini({
    # Potravitug disposition. 'Fit to PK data' / 'Fit to clinical data' in Table S1.
    lvc <- log(3.19); label("Central (blood) volume Vc (L)") # Table S1 row Vc = 3.19 L, fit to PK data
    lvp <- fixed(log(13)); label("Peripheral (other tissue) volume Vp (L)") # Table S1 row Vp = 13 L (Shah and Betts 2012)
    lv_kidney <- fixed(log(0.005)); label("Kidney interstitial space volume VKIS (L)") # Table S1 row VKIS = 0.005 L (Shah and Betts 2012)
    lthalf_ab <- log(31.4); label("Potravitug catabolic half-life, all compartments (day)") # Table S1 row Potravitug t1/2 = 31.4 d, fit to clinical data
    ltdist_ab <- log(52.4); label("Potravitug inter-compartmental distribution half-life Tdist (h)") # Table S1 row Tdist Potravitug = 52.4 (unit not printed; hours, see vignette)
    lpdist_ab <- log(0.441); label("Potravitug tissue:blood partition coefficient Pdist (unitless)") # Table S1 row Pdist Potravitug = 0.441, fit to clinical data

    # VP1 (BKPyV major capsid protein) baseline and turnover
    vp1_per_virion <- fixed(360); label("VP1 copies per BKPyV particle (count)") # Table S1 row VP1 per virus = 360 (Stehle 1996)
    titer_blood <- fixed(7e6); label("Baseline BKPyV titer in blood (virions/L)") # Table S1 row BKPyV blood titer = 7.00E+06 1/L (Hariharan 2006)
    titer_kidney <- fixed(7e9); label("Baseline BKPyV titer in kidney interstitial space (virions/L)") # Table S1 row BKPyV KIS titer = 7.00E+09 1/L (Theodoropoulos 2013)
    lthalf_vp1 <- fixed(log(1440)); label("VP1 half-life (min)") # Table S1 row VP1 t1/2 = 1440 min (Funk 2008); Model Assumptions t_VP1_halflife = 24 h

    # Potravitug-VP1 binding (1:1)
    kon <- fixed(1e-3); label("Potravitug-VP1 association rate constant (1/nM/s)") # Table S1 row konAb-VP1 = 1.00E-03 1/nM/s (Schlosshauer 2004)
    kd <- fixed(60); label("Potravitug-VP1 equilibrium dissociation constant (pM)") # Table S1 row KdAb-VP1 = 60 pM (internal data)

    # Not reported in the paper: molar mass used to convert potravitug mass to
    # nM for binding. Typical IgG1 value; see vignette Assumptions.
    mw <- fixed(150000); label("Potravitug molar mass (g/mol)") # Non-paper value: typical IgG1 molar mass, ~150 kDa
  })

  model({
    # Avogadro constant (1/mol) converts virion titers to molar VP1.
    n_avogadro <- 6.02214076e23

    vc <- exp(lvc)
    vp <- exp(lvp)
    v_kidney <- exp(lv_kidney)

    # First-order catabolism of free potravitug in every compartment (1/h).
    kel <- log(2) / (exp(lthalf_ab) * 24)

    # Pdist/Tdist convention (Table S1 Model Assumptions): for compartments
    # i (central) and j (tissue), k_ij = ln2/Tdist * Pdist / (Pdist + Vi/Vj)
    # and k_ji = ln2/Tdist / (1 + Pdist * Vj/Vi). These are amount-based
    # rate constants; at equilibrium C_j / C_i = Pdist and k_ij + k_ji =
    # ln2/Tdist. The same Pdist and Tdist apply to peripheral and kidney.
    tdist <- exp(ltdist_ab)
    pdist <- exp(lpdist_ab)
    k_cp <- log(2) / tdist * pdist / (pdist + vc / vp)
    k_pc <- log(2) / tdist / (1 + pdist * vp / vc)
    k_ck <- log(2) / tdist * pdist / (pdist + vc / v_kidney)
    k_kc <- log(2) / tdist / (1 + pdist * v_kidney / vc)

    # VP1 turnover. Baseline VP1 (nM) = titer (1/L) * VP1 per virion / N_A
    # * 1e9. Table S1 states that kidney-to-blood VP1 synthesis and
    # distribution use an analytical solution holding the kidney VP1
    # concentration 1000-fold above blood; that net supply is applied here
    # as a zero-order input into each compartment equal to its baseline
    # degradation (see vignette Assumptions).
    kdeg_vp1 <- log(2) / (exp(lthalf_vp1) / 60)
    vp1_blood0 <- titer_blood * vp1_per_virion / n_avogadro * 1e9
    vp1_kidney0 <- titer_kidney * vp1_per_virion / n_avogadro * 1e9
    ksyn_blood <- kdeg_vp1 * vp1_blood0 * vc
    ksyn_kidney <- kdeg_vp1 * vp1_kidney0 * v_kidney

    # Binding: kon from 1/nM/s to 1/nM/h; koff = Kd * kon with Kd in nM.
    kon_h <- kon * 3600
    koff <- kd / 1000 * kon_h

    # Free potravitug molar concentrations (nM); mg/L * 1e6 / (g/mol) = nM.
    ab_c_nm <- central / vc * 1e6 / mw
    ab_k_nm <- kidney / v_kidney * 1e6 / mw

    # Binding fluxes (nmol/h)
    bind_c <- kon_h * ab_c_nm * target - koff * complex
    bind_k <- kon_h * ab_k_nm * target_kidney - koff * complex_kidney

    # Potravitug amounts (mg); bound drug leaves the free pool (nmol -> mg).
    d/dt(central) <- -kel * central - k_cp * central + k_pc * peripheral1 -
      k_ck * central + k_kc * kidney - bind_c * mw / 1e6
    d/dt(peripheral1) <- k_cp * central - k_pc * peripheral1 - kel * peripheral1
    d/dt(kidney) <- k_ck * central - k_kc * kidney - kel * kidney - bind_k * mw / 1e6

    # VP1 and complex amounts (nmol); complex cleared at the VP1 rate.
    d/dt(target) <- ksyn_blood - kdeg_vp1 * target - bind_c
    d/dt(complex) <- bind_c - kdeg_vp1 * complex
    d/dt(target_kidney) <- ksyn_kidney - kdeg_vp1 * target_kidney - bind_k
    d/dt(complex_kidney) <- bind_k - kdeg_vp1 * complex_kidney

    target(0) <- vp1_blood0 * vc
    target_kidney(0) <- vp1_kidney0 * v_kidney

    # Outputs: free potravitug in blood and kidney interstitial space
    # (ug/mL), and VP1 receptor occupancy (% of total VP1 bound).
    Cc <- central / vc
    Ckidney <- kidney / v_kidney
    ro_blood <- 100 * complex / (complex + target)
    ro_kidney <- 100 * complex_kidney / (complex_kidney + target_kidney)
  })
}
