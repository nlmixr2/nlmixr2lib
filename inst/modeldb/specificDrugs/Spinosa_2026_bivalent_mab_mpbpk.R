Spinosa_2026_bivalent_mab_mpbpk <- function() {
  description <- paste(
    "QSP / minimal-PBPK. Generalized bispecific-antibody minimal PBPK model",
    "(Spinosa 2026, Case Study 4) in its bivalent monospecific configuration:",
    "a prototypical IgG that binds the therapeutic receptor with both arms",
    "(10 nM), the comparator for the tumor-targeted bispecific antibody in",
    "patients with solid tumors. Second-generation Cao-Jusko minimal PBPK",
    "backbone (plasma, leaky and tight interstitial tissue, lymph) plus a tumor",
    "interstitial compartment fed by 7% of lymph flow; nonspecific clearance",
    "from plasma. The therapeutic receptor is expressed at the same level",
    "(10^4 per cell, 10^6 cells/mL) in plasma, leaky, tight and tumor",
    "compartments. First binding uses twice the monovalent on-rate (two",
    "identical arms); the single-bound antibody cross-links a second receptor",
    "on the same cell with an effective avidity factor on the on-rate.",
    "Receptor-bound antibody is internalized at the receptor turnover rate and",
    "a doubly bound antibody internalizes through either receptor, recycling",
    "the other. A hypothetical molecule: all drug parameters are design values."
  )
  reference <- paste(
    "Spinosa P, Joslyn L, Ramanujan S, Gadkar K, Hosseini I. A Generalized",
    "Minimal PBPK-PD Model of Bispecific Antibodies: Case Studies and",
    "Applications in Drug Development. CPT Pharmacometrics Syst Pharmacol.",
    "2026;15(2):e70167. doi:10.1002/psp4.70167. Model equations from the",
    "Supporting Information Data S1 (SimBiology project PSP-2025-0197-s04.sbproj,",
    "model 'm1', FLAG_BIVALENT_BINDING = 1); biology from the parameter workbook",
    "PSP-2025-0197-s07.xlsx sheet 'Case Study 4' and Table S3; antibody",
    "parameters from Table S4. The tumor-targeted bispecific it is compared",
    "with is modellib('Spinosa_2026_tumor_targeted_bsab_mpbpk').",
    sep = " "
  )
  vignette <- "Spinosa_2026_bispecific_mpbpk"

  # Every state holds a CONCENTRATION (nM) in a body-weight-normalised volume
  # (mL/kg), as in the published SimBiology model. Doses are IV in mg/kg and
  # f(central) converts mg/kg into nM (g/L divided by g/mol, times 1e9).
  dosing <- c("central")
  units <- list(time = "day", dosing = "mg/kg", concentration = "ug/mL")

  # drug_ther_<location> is antibody bound by one arm, drug2_ther_<location>
  # antibody bound by both arms to two receptors on the same cell.
  paper_specific_compartments <- c(
    "tumor_isf",
    "target_ther_central",
    "target_ther_leaky",
    "target_ther_tight",
    "target_ther_tumor",
    "drug_ther_central",
    "drug_ther_leaky",
    "drug_ther_tight",
    "drug_ther_tumor",
    "drug2_ther_central",
    "drug2_ther_leaky",
    "drug2_ther_tight",
    "drug2_ther_tumor"
  )

  compartmentData <- list(
    central = list(analyte = "bivalent antibody (free)", units = "nM", specimen = "plasma", verified = TRUE),
    leaky = list(analyte = "bivalent antibody (free)", units = "nM", specimen = "tissue", verified = TRUE),
    tight = list(analyte = "bivalent antibody (free)", units = "nM", specimen = "tissue", verified = TRUE),
    tumor_isf = list(analyte = "bivalent antibody (free)", units = "nM", specimen = "tumor", verified = TRUE),
    lymph = list(analyte = "bivalent antibody (free)", units = "nM", specimen = "lymph", verified = TRUE),
    target_ther_central = list(
      analyte = "therapeutic receptor (free)",
      units = "nM",
      specimen = "whole blood",
      verified = TRUE
    ),
    target_ther_leaky = list(
      analyte = "therapeutic receptor (free)",
      units = "nM",
      specimen = "tissue",
      verified = TRUE
    ),
    target_ther_tight = list(
      analyte = "therapeutic receptor (free)",
      units = "nM",
      specimen = "tissue",
      verified = TRUE
    ),
    target_ther_tumor = list(
      analyte = "therapeutic receptor (free)",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    drug_ther_central = list(
      analyte = "antibody bound by one arm",
      units = "nM",
      specimen = "whole blood",
      verified = TRUE
    ),
    drug_ther_leaky = list(analyte = "antibody bound by one arm", units = "nM", specimen = "tissue", verified = TRUE),
    drug_ther_tight = list(analyte = "antibody bound by one arm", units = "nM", specimen = "tissue", verified = TRUE),
    drug_ther_tumor = list(analyte = "antibody bound by one arm", units = "nM", specimen = "tumor", verified = TRUE),
    drug2_ther_central = list(
      analyte = "antibody bound by both arms",
      units = "nM",
      specimen = "whole blood",
      verified = TRUE
    ),
    drug2_ther_leaky = list(
      analyte = "antibody bound by both arms",
      units = "nM",
      specimen = "tissue",
      verified = TRUE
    ),
    drug2_ther_tight = list(
      analyte = "antibody bound by both arms",
      units = "nM",
      specimen = "tissue",
      verified = TRUE
    ),
    drug2_ther_tumor = list(analyte = "antibody bound by both arms", units = "nM", specimen = "tumor", verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species = "human",
    n_subjects = NA_integer_,
    n_studies = 0L,
    age_range = "not applicable (simulation case study)",
    weight_range = "70 kg reference body weight (parameter workbook, Case Study 4); not used for mg/kg dosing",
    sex_female_pct = NA_real_,
    race_ethnicity = "not applicable",
    disease_state = "patients with solid tumors (simulated; no clinical data)",
    dose_range = "0.01-10 mg/kg single IV dose (Figure 7e,f)",
    regions = "not applicable",
    notes = paste(
      "A design-exploration comparator, not a fit: prototypical receptor biology",
      "(Table S3) and antibody parameters (Table S4)."
    )
  )

  ini({
    # Antibody disposition (Table S4; workbook sheet 'Case Study 4')
    lcl <- fixed(log(2.8272)); label("Nonspecific clearance from plasma (mL/day/kg)") # workbook CLp_D1 = 2.8272; Table S4 = 2.8 mL/kg/day
    lvc <- fixed(log(37.14)); label("Plasma (central) volume (mL/kg)") # workbook V_cen = 37.14

    # Minimal PBPK physiology with a tumor compartment (workbook)
    lvleaky <- fixed(log(62.4)); label("Leaky-tissue interstitial volume (mL/kg)") # workbook V_lea = 62.4
    lvtight <- fixed(log(115.9)); label("Tight-tissue interstitial volume (mL/kg)") # workbook V_tig = 115.9
    lvtumor <- fixed(log(4)); label("Tumor interstitial volume (mL/kg)") # workbook V_eff = 4
    lvlymph <- fixed(log(37.14)); label("Lymph volume (mL/kg)") # workbook V_lym = 37.14
    llymphflow <- fixed(log(41.4288)); label("Total lymph flow (mL/day/kg)") # workbook L_D1 = 41.4288
    flymph_leaky <- fixed(0.58); label("Fraction of lymph flow through leaky tissue (fraction)") # workbook f_lea_D1 = 0.58
    flymph_tight <- fixed(0.35); label("Fraction of lymph flow through tight tissue (fraction)") # workbook f_tig_D1 = 0.35
    flymph_tumor <- fixed(0.07); label("Fraction of lymph flow through the tumor (fraction)") # workbook f_eff_D1 = 0.07
    sigma_leaky <- fixed(0.7); label("Vascular reflection coefficient, leaky tissue (unitless)") # workbook sig_lea_D1 = 0.7
    sigma_tight <- fixed(0.95); label("Vascular reflection coefficient, tight tissue (unitless)") # workbook sig_tig_D1 = 0.95
    sigma_tumor <- fixed(0.7); label("Vascular reflection coefficient, tumor (unitless)") # workbook sig_eff_D1 = 0.7
    sigma_l <- fixed(0.2); label("Lymphatic reflection coefficient (unitless)") # workbook sig_lym_D1 = 0.2
    mw <- fixed(150000); label("Antibody molecular weight (g/mol)") # workbook MWab_D1 = 150000

    # Binding
    kon_ther <- fixed(3.6); label("Monovalent association rate constant for the therapeutic receptor (1/nM/day)") # workbook kon_R1_D1 = 3.6
    kd_ther <- fixed(10); label("Equilibrium dissociation constant for the therapeutic receptor (nM)") # Table S4 bivalent antibody = 10 nM
    chi_ther <- fixed(1000); label("Effective avidity factor for the second arm binding (unitless)") # project variant '6. AVIDITY: same cell bv mAb' avid_R1 = 1e9 per 1e6 cells/mL

    # Receptor expression (Table S3) and turnover (workbook)
    cells_per_ml <- fixed(1e6); label("Cells expressing the therapeutic receptor in every compartment (cells/mL)") # Table S3 = 10^6 cells/mL
    rpc_ther <- fixed(1e4); label("Therapeutic receptors per cell, all compartments (receptors/cell)") # Table S3 = 10^4; workbook rpc_mR1_* = 10000
    thalf_ther <- fixed(1); label("Therapeutic receptor turnover half-life (day)") # workbook Thalf_mR1 = 1

    propSd <- fixed(0); label("Proportional residual error (fraction; not estimated in the source)")
  })

  model({
    n_avogadro <- 6.022e23
    cl <- exp(lcl)
    vc <- exp(lvc)
    vleaky <- exp(lvleaky)
    vtight <- exp(lvtight)
    vtumor <- exp(lvtumor)
    vlymph <- exp(lvlymph)
    lymphflow <- exp(llymphflow)
    l_leaky <- flymph_leaky * lymphflow
    l_tight <- flymph_tight * lymphflow
    l_tumor <- flymph_tumor * lymphflow

    k_c_leaky <- l_leaky * (1 - sigma_leaky) / vc
    k_c_tight <- l_tight * (1 - sigma_tight) / vc
    k_c_tumor <- l_tumor * (1 - sigma_tumor) / vc
    k_leaky_l <- l_leaky * (1 - sigma_l) / vleaky
    k_tight_l <- l_tight * (1 - sigma_l) / vtight
    k_tumor_l <- l_tumor * (1 - sigma_l) / vtumor
    k_l_c <- lymphflow / vlymph
    kel <- cl / vc

    bl_ther <- cells_per_ml * rpc_ther / n_avogadro * 1e12
    kdeg_ther <- log(2) / thalf_ther
    # Two identical arms: first binding on-rate doubled (published
    # kon_R1_homodimer = 2 * kon_R1); the doubly bound antibody dissociates
    # through either arm, so its off-rate is doubled as well.
    kon2_ther <- 2 * kon_ther

    target_ther_central(0) <- bl_ther
    target_ther_leaky(0) <- bl_ther
    target_ther_tight(0) <- bl_ther
    target_ther_tumor(0) <- bl_ther

    b1_c <- kon2_ther * central * target_ther_central - kon_ther * kd_ther * drug_ther_central
    b2_c <- chi_ther * kon_ther * drug_ther_central * target_ther_central - kon2_ther * kd_ther * drug2_ther_central
    b1_lk <- kon2_ther * leaky * target_ther_leaky - kon_ther * kd_ther * drug_ther_leaky
    b2_lk <- chi_ther * kon_ther * drug_ther_leaky * target_ther_leaky - kon2_ther * kd_ther * drug2_ther_leaky
    b1_tt <- kon2_ther * tight * target_ther_tight - kon_ther * kd_ther * drug_ther_tight
    b2_tt <- chi_ther * kon_ther * drug_ther_tight * target_ther_tight - kon2_ther * kd_ther * drug2_ther_tight
    b1_tu <- kon2_ther * tumor_isf * target_ther_tumor - kon_ther * kd_ther * drug_ther_tumor
    b2_tu <- chi_ther * kon_ther * drug_ther_tumor * target_ther_tumor - kon2_ther * kd_ther * drug2_ther_tumor

    d/dt(central) <- -kel * central - (k_c_leaky + k_c_tight + k_c_tumor) * central +
      k_l_c * vlymph / vc * lymph - b1_c
    d/dt(leaky) <- k_c_leaky * vc / vleaky * central - k_leaky_l * leaky - b1_lk
    d/dt(tight) <- k_c_tight * vc / vtight * central - k_tight_l * tight - b1_tt
    d/dt(tumor_isf) <- k_c_tumor * vc / vtumor * central - k_tumor_l * tumor_isf - b1_tu
    d/dt(lymph) <- k_leaky_l * vleaky / vlymph * leaky + k_tight_l * vtight / vlymph * tight +
      k_tumor_l * vtumor / vlymph * tumor_isf - k_l_c * lymph

    # A doubly bound antibody internalizes through either of its two receptors
    # (rate 2 * kdeg) and returns the other receptor to the free pool.
    d/dt(target_ther_central) <- kdeg_ther * bl_ther - kdeg_ther * target_ther_central - b1_c - b2_c + 2 * kdeg_ther * drug2_ther_central
    d/dt(target_ther_leaky) <- kdeg_ther * bl_ther - kdeg_ther * target_ther_leaky - b1_lk - b2_lk + 2 * kdeg_ther * drug2_ther_leaky
    d/dt(target_ther_tight) <- kdeg_ther * bl_ther - kdeg_ther * target_ther_tight - b1_tt - b2_tt + 2 * kdeg_ther * drug2_ther_tight
    d/dt(target_ther_tumor) <- kdeg_ther * bl_ther - kdeg_ther * target_ther_tumor - b1_tu - b2_tu + 2 * kdeg_ther * drug2_ther_tumor

    d/dt(drug_ther_central) <- b1_c - b2_c - kdeg_ther * drug_ther_central
    d/dt(drug_ther_leaky) <- b1_lk - b2_lk - kdeg_ther * drug_ther_leaky
    d/dt(drug_ther_tight) <- b1_tt - b2_tt - kdeg_ther * drug_ther_tight
    d/dt(drug_ther_tumor) <- b1_tu - b2_tu - kdeg_ther * drug_ther_tumor

    d/dt(drug2_ther_central) <- b2_c - 2 * kdeg_ther * drug2_ther_central
    d/dt(drug2_ther_leaky) <- b2_lk - 2 * kdeg_ther * drug2_ther_leaky
    d/dt(drug2_ther_tight) <- b2_tt - 2 * kdeg_ther * drug2_ther_tight
    d/dt(drug2_ther_tumor) <- b2_tu - 2 * kdeg_ther * drug2_ther_tumor

    f(central) <- 1e9 / (mw * vc)

    Cc <- central * mw / 1e6

    ro_ther_central <- 100 * (1 - target_ther_central / bl_ther)
    ro_ther_leaky <- 100 * (1 - target_ther_leaky / bl_ther)
    ro_ther_tight <- 100 * (1 - target_ther_tight / bl_ther)
    ro_ther_tumor <- 100 * (1 - target_ther_tumor / bl_ther)

    Cc ~ prop(propSd)
  })
}
