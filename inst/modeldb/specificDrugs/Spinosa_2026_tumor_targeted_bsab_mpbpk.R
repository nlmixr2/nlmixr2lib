Spinosa_2026_tumor_targeted_bsab_mpbpk <- function() {
  description <- paste(
    "QSP / minimal-PBPK. Generalized bispecific-antibody minimal PBPK model",
    "(Spinosa 2026, Case Study 4) parameterized for a prototypical",
    "tumor-targeted bispecific antibody in patients with solid tumors: one arm",
    "binds a therapeutic receptor with low affinity (100 nM) and the other a",
    "tumor-enriched targeting receptor with high affinity (0.05 nM). Second-",
    "generation Cao-Jusko minimal PBPK backbone (plasma, leaky and tight",
    "interstitial tissue, lymph) plus a tumor (efficacy) interstitial",
    "compartment fed by 7% of lymph flow; nonspecific clearance from plasma.",
    "Both receptors are expressed on the same cells (10^6 cells/mL) in plasma,",
    "leaky, tight and tumor compartments; the targeting receptor is 100-fold",
    "enriched in the tumor. The drug binds either receptor and the single-bound",
    "complex binds the other receptor on the same cell (cis-binding) with an",
    "effective avidity factor multiplying the second on-rate. Receptor-bound",
    "drug is internalized at the receptor turnover rate. The model output of",
    "interest is therapeutic-receptor occupancy in blood (safety) versus tumor",
    "(efficacy). A hypothetical molecule: all drug parameters are design values."
  )
  reference <- paste(
    "Spinosa P, Joslyn L, Ramanujan S, Gadkar K, Hosseini I. A Generalized",
    "Minimal PBPK-PD Model of Bispecific Antibodies: Case Studies and",
    "Applications in Drug Development. CPT Pharmacometrics Syst Pharmacol.",
    "2026;15(2):e70167. doi:10.1002/psp4.70167. Model equations and default",
    "parameters from the Supporting Information Data S1 (SimBiology project",
    "PSP-2025-0197-s04.sbproj, model 'm1'; parameter workbook",
    "PSP-2025-0197-s07.xlsx sheet 'Case Study 4'; Tables S3 and S4).",
    "The bivalent monospecific comparator is modellib('Spinosa_2026_bivalent_mab_mpbpk').",
    sep = " "
  )
  vignette <- "Spinosa_2026_bispecific_mpbpk"

  # Every state holds a CONCENTRATION (nM) in a body-weight-normalised volume
  # (mL/kg), as in the published SimBiology model. Doses are IV in mg/kg and
  # f(central) converts mg/kg into nM (g/L divided by g/mol, times 1e9).
  dosing <- c("central")
  units <- list(time = "day", dosing = "mg/kg", concentration = "ug/mL")

  # `tumor_isf` is the free drug in the tumor (the model's 'efficacy'
  # compartment, an interstitial space like `leaky` and `tight`). Receptor
  # species are named <species>_<receptor>_<location>, with `ther` the
  # therapeutic receptor and `tgt` the tumor-targeting receptor, and
  # trimer_<location> the drug bound to both receptors on one cell.
  paper_specific_compartments <- c(
    "tumor_isf",
    "target_ther_central",
    "target_ther_leaky",
    "target_ther_tight",
    "target_ther_tumor",
    "target_tgt_central",
    "target_tgt_leaky",
    "target_tgt_tight",
    "target_tgt_tumor",
    "drug_ther_central",
    "drug_ther_leaky",
    "drug_ther_tight",
    "drug_ther_tumor",
    "drug_tgt_central",
    "drug_tgt_leaky",
    "drug_tgt_tight",
    "drug_tgt_tumor",
    "trimer_central",
    "trimer_leaky",
    "trimer_tight",
    "trimer_tumor"
  )

  compartmentData <- list(
    central = list(analyte = "bispecific antibody (free)", units = "nM", specimen = "plasma", verified = TRUE),
    leaky = list(analyte = "bispecific antibody (free)", units = "nM", specimen = "tissue", verified = TRUE),
    tight = list(analyte = "bispecific antibody (free)", units = "nM", specimen = "tissue", verified = TRUE),
    tumor_isf = list(analyte = "bispecific antibody (free)", units = "nM", specimen = "tumor", verified = TRUE),
    lymph = list(analyte = "bispecific antibody (free)", units = "nM", specimen = "lymph", verified = TRUE),
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
    target_tgt_central = list(
      analyte = "targeting receptor (free)",
      units = "nM",
      specimen = "whole blood",
      verified = TRUE
    ),
    target_tgt_leaky = list(analyte = "targeting receptor (free)", units = "nM", specimen = "tissue", verified = TRUE),
    target_tgt_tight = list(analyte = "targeting receptor (free)", units = "nM", specimen = "tissue", verified = TRUE),
    target_tgt_tumor = list(analyte = "targeting receptor (free)", units = "nM", specimen = "tumor", verified = TRUE),
    drug_ther_central = list(
      analyte = "drug-therapeutic receptor dimer",
      units = "nM",
      specimen = "whole blood",
      verified = TRUE
    ),
    drug_ther_leaky = list(
      analyte = "drug-therapeutic receptor dimer",
      units = "nM",
      specimen = "tissue",
      verified = TRUE
    ),
    drug_ther_tight = list(
      analyte = "drug-therapeutic receptor dimer",
      units = "nM",
      specimen = "tissue",
      verified = TRUE
    ),
    drug_ther_tumor = list(
      analyte = "drug-therapeutic receptor dimer",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    ),
    drug_tgt_central = list(
      analyte = "drug-targeting receptor dimer",
      units = "nM",
      specimen = "whole blood",
      verified = TRUE
    ),
    drug_tgt_leaky = list(
      analyte = "drug-targeting receptor dimer",
      units = "nM",
      specimen = "tissue",
      verified = TRUE
    ),
    drug_tgt_tight = list(
      analyte = "drug-targeting receptor dimer",
      units = "nM",
      specimen = "tissue",
      verified = TRUE
    ),
    drug_tgt_tumor = list(analyte = "drug-targeting receptor dimer", units = "nM", specimen = "tumor", verified = TRUE),
    trimer_central = list(
      analyte = "therapeutic receptor-drug-targeting receptor trimer",
      units = "nM",
      specimen = "whole blood",
      verified = TRUE
    ),
    trimer_leaky = list(
      analyte = "therapeutic receptor-drug-targeting receptor trimer",
      units = "nM",
      specimen = "tissue",
      verified = TRUE
    ),
    trimer_tight = list(
      analyte = "therapeutic receptor-drug-targeting receptor trimer",
      units = "nM",
      specimen = "tissue",
      verified = TRUE
    ),
    trimer_tumor = list(
      analyte = "therapeutic receptor-drug-targeting receptor trimer",
      units = "nM",
      specimen = "tumor",
      verified = TRUE
    )
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
    dose_range = "0.01-10 mg/kg single IV dose (Figure 7)",
    regions = "not applicable",
    notes = paste(
      "A design-exploration case study, not a fit: the receptor expression levels",
      "(Table S3), affinities and clearance (Table S4) are prototypical values chosen",
      "by the authors to study how the affinity ratio of the two arms and the",
      "effective avidity shape blood-versus-tumor therapeutic-receptor occupancy."
    )
  )

  ini({
    # Drug disposition (Table S4; workbook sheet 'Case Study 4')
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

    # Binding (R1 = therapeutic receptor, R2 = targeting receptor)
    kon_ther <- fixed(3.6); label("Association rate constant for the therapeutic receptor (1/nM/day)") # workbook kon_R1_D1 = 3.6
    kon_tgt <- fixed(3.6); label("Association rate constant for the targeting receptor (1/nM/day)") # workbook kon_R2_D1 = 3.6
    kd_ther <- fixed(100); label("Equilibrium dissociation constant for the therapeutic receptor (nM)") # Table S4 = 100 nM; workbook Kd_R1_D1 = 100
    kd_tgt <- fixed(0.05); label("Equilibrium dissociation constant for the targeting receptor (nM)") # Table S4 = 0.05 nM; workbook Kd_R2_D1 = 0.05
    chi_ther <- fixed(1000); label("Effective avidity factor for the therapeutic arm binding second (unitless)") # Figure 7 caption: avidity factor held at 1000; workbook avid_R1_D1 = 1e9 per 1e6 cells/mL
    chi_tgt <- fixed(1e-6); label("Effective avidity factor for the targeting arm binding second (unitless)") # workbook avid_R2_D1 = 1 per 1e6 cells/mL

    # Receptor expression (Table S3) and turnover (workbook)
    cells_per_ml <- fixed(1e6); label("Cells expressing both receptors in every compartment (cells/mL)") # Table S3 = 10^6 cells/mL; workbook mR1/mR2_cellspermL = 1e6
    rpc_ther <- fixed(1e4); label("Therapeutic receptors per cell, all compartments (receptors/cell)") # Table S3 = 10^4; workbook rpc_mR1_* = 10000
    rpc_tgt <- fixed(1e3); label("Targeting receptors per cell, non-tumor compartments (receptors/cell)") # Table S3 = 10^3; workbook rpc_mR2_cen/lea/tig = 1000
    rpc_tgt_tumor <- fixed(1e5); label("Targeting receptors per cell in the tumor (receptors/cell)") # Table S3 = 10^5; workbook rpc_mR2_eff = 1e5
    thalf_ther <- fixed(1); label("Therapeutic receptor turnover half-life (day)") # workbook Thalf_mR1 = 1
    thalf_tgt <- fixed(1); label("Targeting receptor turnover half-life (day)") # workbook Thalf_mR2 = 1

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

    # Receptors per mL -> nM: cells/mL * receptors/cell / Avogadro * 1e3 mL/L * 1e9 nmol/mol
    bl_ther <- cells_per_ml * rpc_ther / n_avogadro * 1e12
    bl_tgt <- cells_per_ml * rpc_tgt / n_avogadro * 1e12
    bl_tgt_tumor <- cells_per_ml * rpc_tgt_tumor / n_avogadro * 1e12

    # Receptor-bound drug internalizes at the receptor turnover rate (fold
    # change in internalization = 1 in the workbook).
    kdeg_ther <- log(2) / thalf_ther
    kdeg_tgt <- log(2) / thalf_tgt

    target_ther_central(0) <- bl_ther
    target_ther_leaky(0) <- bl_ther
    target_ther_tight(0) <- bl_ther
    target_ther_tumor(0) <- bl_ther
    target_tgt_central(0) <- bl_tgt
    target_tgt_leaky(0) <- bl_tgt
    target_tgt_tight(0) <- bl_tgt
    target_tgt_tumor(0) <- bl_tgt_tumor

    # Binding fluxes (nM/day) per location: first binding to either receptor,
    # then cis-binding of the single-bound complex to the other receptor with
    # its on-rate multiplied by the effective avidity factor.
    bth_c <- kon_ther * (central * target_ther_central - kd_ther * drug_ther_central)
    btg_c <- kon_tgt * (central * target_tgt_central - kd_tgt * drug_tgt_central)
    bth2tg_c <- chi_tgt * kon_tgt * drug_ther_central * target_tgt_central - kon_tgt * kd_tgt * trimer_central
    btg2th_c <- chi_ther * kon_ther * drug_tgt_central * target_ther_central - kon_ther * kd_ther * trimer_central

    bth_lk <- kon_ther * (leaky * target_ther_leaky - kd_ther * drug_ther_leaky)
    btg_lk <- kon_tgt * (leaky * target_tgt_leaky - kd_tgt * drug_tgt_leaky)
    bth2tg_lk <- chi_tgt * kon_tgt * drug_ther_leaky * target_tgt_leaky - kon_tgt * kd_tgt * trimer_leaky
    btg2th_lk <- chi_ther * kon_ther * drug_tgt_leaky * target_ther_leaky - kon_ther * kd_ther * trimer_leaky

    bth_tt <- kon_ther * (tight * target_ther_tight - kd_ther * drug_ther_tight)
    btg_tt <- kon_tgt * (tight * target_tgt_tight - kd_tgt * drug_tgt_tight)
    bth2tg_tt <- chi_tgt * kon_tgt * drug_ther_tight * target_tgt_tight - kon_tgt * kd_tgt * trimer_tight
    btg2th_tt <- chi_ther * kon_ther * drug_tgt_tight * target_ther_tight - kon_ther * kd_ther * trimer_tight

    bth_tu <- kon_ther * (tumor_isf * target_ther_tumor - kd_ther * drug_ther_tumor)
    btg_tu <- kon_tgt * (tumor_isf * target_tgt_tumor - kd_tgt * drug_tgt_tumor)
    bth2tg_tu <- chi_tgt * kon_tgt * drug_ther_tumor * target_tgt_tumor - kon_tgt * kd_tgt * trimer_tumor
    btg2th_tu <- chi_ther * kon_ther * drug_tgt_tumor * target_ther_tumor - kon_ther * kd_ther * trimer_tumor

    # Free drug
    d/dt(central) <- -kel * central - (k_c_leaky + k_c_tight + k_c_tumor) * central +
      k_l_c * vlymph / vc * lymph - bth_c - btg_c
    d/dt(leaky) <- k_c_leaky * vc / vleaky * central - k_leaky_l * leaky - bth_lk - btg_lk
    d/dt(tight) <- k_c_tight * vc / vtight * central - k_tight_l * tight - bth_tt - btg_tt
    d/dt(tumor_isf) <- k_c_tumor * vc / vtumor * central - k_tumor_l * tumor_isf - bth_tu - btg_tu
    d/dt(lymph) <- k_leaky_l * vleaky / vlymph * leaky + k_tight_l * vtight / vlymph * tight +
      k_tumor_l * vtumor / vlymph * tumor_isf - k_l_c * lymph

    # Free receptors: synthesis, turnover, binding, and recycling of the
    # partner receptor when a trimer internalizes through the other one.
    d/dt(target_ther_central) <- kdeg_ther * bl_ther - kdeg_ther * target_ther_central - bth_c - btg2th_c + kdeg_tgt * trimer_central
    d/dt(target_ther_leaky) <- kdeg_ther * bl_ther - kdeg_ther * target_ther_leaky - bth_lk - btg2th_lk + kdeg_tgt * trimer_leaky
    d/dt(target_ther_tight) <- kdeg_ther * bl_ther - kdeg_ther * target_ther_tight - bth_tt - btg2th_tt + kdeg_tgt * trimer_tight
    d/dt(target_ther_tumor) <- kdeg_ther * bl_ther - kdeg_ther * target_ther_tumor - bth_tu - btg2th_tu + kdeg_tgt * trimer_tumor
    d/dt(target_tgt_central) <- kdeg_tgt * bl_tgt - kdeg_tgt * target_tgt_central - btg_c - bth2tg_c + kdeg_ther * trimer_central
    d/dt(target_tgt_leaky) <- kdeg_tgt * bl_tgt - kdeg_tgt * target_tgt_leaky - btg_lk - bth2tg_lk + kdeg_ther * trimer_leaky
    d/dt(target_tgt_tight) <- kdeg_tgt * bl_tgt - kdeg_tgt * target_tgt_tight - btg_tt - bth2tg_tt + kdeg_ther * trimer_tight
    d/dt(target_tgt_tumor) <- kdeg_tgt * bl_tgt_tumor - kdeg_tgt * target_tgt_tumor - btg_tu - bth2tg_tu + kdeg_ther * trimer_tumor

    # Single-bound dimers (internalized at the receptor turnover rate)
    d/dt(drug_ther_central) <- bth_c - bth2tg_c - kdeg_ther * drug_ther_central
    d/dt(drug_ther_leaky) <- bth_lk - bth2tg_lk - kdeg_ther * drug_ther_leaky
    d/dt(drug_ther_tight) <- bth_tt - bth2tg_tt - kdeg_ther * drug_ther_tight
    d/dt(drug_ther_tumor) <- bth_tu - bth2tg_tu - kdeg_ther * drug_ther_tumor
    d/dt(drug_tgt_central) <- btg_c - btg2th_c - kdeg_tgt * drug_tgt_central
    d/dt(drug_tgt_leaky) <- btg_lk - btg2th_lk - kdeg_tgt * drug_tgt_leaky
    d/dt(drug_tgt_tight) <- btg_tt - btg2th_tt - kdeg_tgt * drug_tgt_tight
    d/dt(drug_tgt_tumor) <- btg_tu - btg2th_tu - kdeg_tgt * drug_tgt_tumor

    # Trimers internalize through each receptor at that receptor's rate
    d/dt(trimer_central) <- bth2tg_c + btg2th_c - (kdeg_ther + kdeg_tgt) * trimer_central
    d/dt(trimer_leaky) <- bth2tg_lk + btg2th_lk - (kdeg_ther + kdeg_tgt) * trimer_leaky
    d/dt(trimer_tight) <- bth2tg_tt + btg2th_tt - (kdeg_ther + kdeg_tgt) * trimer_tight
    d/dt(trimer_tumor) <- bth2tg_tu + btg2th_tu - (kdeg_ther + kdeg_tgt) * trimer_tumor

    f(central) <- 1e9 / (mw * vc)

    # Free drug in plasma
    Cc <- central * mw / 1e6

    # Therapeutic-receptor occupancy, 100 * (1 - free / baseline), as in the
    # published model's RO_mR1_<location>
    ro_ther_central <- 100 * (1 - target_ther_central / bl_ther)
    ro_ther_leaky <- 100 * (1 - target_ther_leaky / bl_ther)
    ro_ther_tight <- 100 * (1 - target_ther_tight / bl_ther)
    ro_ther_tumor <- 100 * (1 - target_ther_tumor / bl_ther)

    Cc ~ prop(propSd)
  })
}
