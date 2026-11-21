Spinosa_2026_mosunetuzumab_cyno_mpbpk <- function() {
  description <- paste(
    "Preclinical (cynomolgus monkey). QSP / minimal-PBPK. Generalized",
    "bispecific-antibody minimal PBPK model (Spinosa 2026, Case Study 2)",
    "parameterized for mosunetuzumab, a CD20 x CD3 T-cell-dependent bispecific",
    "antibody, in cynomolgus monkeys. Second-generation Cao-Jusko minimal PBPK",
    "backbone (plasma, leaky and tight interstitial tissue, lymph) with",
    "nonspecific clearance from plasma. CD20 (B cells) and CD3 (T cells) are",
    "expressed in the central compartment only, each with zero-order synthesis",
    "and first-order turnover. The drug binds either receptor and the",
    "single-bound complex binds the opposite receptor on the other cell",
    "(trans-binding) to form the drug-CD20-CD3 trimer; receptor-bound drug is",
    "internalized at the receptor turnover rate, and a trimer internalizes",
    "through each receptor at that receptor's rate while the other receptor is",
    "recycled. CD3- and CD20-mediated disposition are therefore explicit target-",
    "mediated clearance routes alongside the nonspecific clearance; cumulative",
    "elimination through each route is tracked to apportion total clearance."
  )
  reference <- paste(
    "Spinosa P, Joslyn L, Ramanujan S, Gadkar K, Hosseini I. A Generalized",
    "Minimal PBPK-PD Model of Bispecific Antibodies: Case Studies and",
    "Applications in Drug Development. CPT Pharmacometrics Syst Pharmacol.",
    "2026;15(2):e70167. doi:10.1002/psp4.70167. Model equations and default",
    "parameters from the Supporting Information Data S1 (SimBiology project",
    "PSP-2025-0197-s04.sbproj, model 'm1'; parameter workbook",
    "PSP-2025-0197-s07.xlsx sheets 'Case Study 2 Mosunetuzumab' and 'Case",
    "Study 3'; Table S2). Cynomolgus monkey mosunetuzumab PK from Hosseini I",
    "et al. NPJ Syst Biol Appl. 2020;6:28 doi:10.1038/s41540-020-00145-7.",
    "Population-PK models of the same drug: modellib('Ferl_2018_mosunetuzumab_cyno').",
    sep = " "
  )
  vignette <- "Spinosa_2026_bispecific_mpbpk"

  # Every state holds a CONCENTRATION (nM) in a body-weight-normalised volume
  # (mL/kg), as in the published SimBiology model. Doses are IV in mg/kg;
  # f(central) converts mg/kg into nM: (mg/kg) / (mL/kg) is g/L, divided by
  # the molecular weight (g/mol) and multiplied by 1e9 gives nM.
  dosing <- c("central")
  units <- list(time = "day", dosing = "mg/kg", concentration = "ug/mL")

  # target_cd3_central / drug_cd3_central / trimer are registered canonicals;
  # the CD20 arm mirrors them. The elim_* states integrate the elimination
  # flux through each route (nmol/kg) so that the fraction of total clearance
  # by route (Figure 6b,c) can be read off a simulation.
  paper_specific_compartments <- c(
    "target_cd20_central",
    "drug_cd20_central",
    "elim_ns",
    "elim_cd20",
    "elim_cd3"
  )

  compartmentData <- list(
    central = list(analyte = "mosunetuzumab (free)", units = "nM", specimen = "serum", verified = TRUE),
    leaky = list(analyte = "mosunetuzumab (free)", units = "nM", specimen = "tissue", verified = TRUE),
    tight = list(analyte = "mosunetuzumab (free)", units = "nM", specimen = "tissue", verified = TRUE),
    lymph = list(analyte = "mosunetuzumab (free)", units = "nM", specimen = "lymph", verified = TRUE),
    target_cd20_central = list(analyte = "CD20 (free)", units = "nM", specimen = "whole blood", verified = TRUE),
    target_cd3_central = list(analyte = "CD3 (free)", units = "nM", specimen = "whole blood", verified = TRUE),
    drug_cd20_central = list(
      analyte = "mosunetuzumab-CD20 dimer",
      units = "nM",
      specimen = "whole blood",
      verified = TRUE
    ),
    drug_cd3_central = list(
      analyte = "mosunetuzumab-CD3 dimer",
      units = "nM",
      specimen = "whole blood",
      verified = TRUE
    ),
    trimer = list(analyte = "CD20-mosunetuzumab-CD3 trimer", units = "nM", specimen = "whole blood", verified = TRUE),
    elim_ns = list(
      analyte = "mosunetuzumab eliminated nonspecifically",
      units = "nmol/kg",
      specimen = "not applicable",
      verified = TRUE
    ),
    elim_cd20 = list(
      analyte = "mosunetuzumab eliminated via CD20",
      units = "nmol/kg",
      specimen = "not applicable",
      verified = TRUE
    ),
    elim_cd3 = list(
      analyte = "mosunetuzumab eliminated via CD3",
      units = "nmol/kg",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  covariateData <- list()

  population <- list(
    species = "cynomolgus monkey",
    n_subjects = NA_integer_,
    n_studies = 1L,
    age_range = "not reported",
    weight_range = "3.12 kg reference body weight (parameter workbook, Case Study 2); not used for mg/kg dosing",
    sex_female_pct = NA_real_,
    race_ethnicity = "not applicable",
    disease_state = "healthy cynomolgus monkeys (nonclinical multi-dose PK study of mosunetuzumab; Hosseini 2020)",
    dose_range = "0.01, 0.1 and 1 mg/kg IV once weekly x 4 (Figure 5d)",
    regions = "not applicable",
    notes = paste(
      "Mechanistic model; only the nonspecific clearance was fitted to the mosunetuzumab",
      "PK. The total CD3 capacity (19 nM) was fitted beforehand to anti-gD/CD3 PK and the",
      "total CD20 capacity (4.4 nM) to obinutuzumab PK in cynomolgus monkeys, each with a",
      "+/-30% range for between-animal heterogeneity (Figure 4 workflow; Table S2).",
      "No inter-individual variability or residual error is estimated."
    )
  )

  ini({
    # Drug disposition (Table S2; workbook sheet 'Case Study 2 Mosunetuzumab')
    lcl <- log(4.0); label("Nonspecific clearance from plasma (mL/day/kg)") # Table S2 CLns = 4.0 mL/kg/day (fitted); workbook CLp_D1 = 4
    lvc <- fixed(log(37.14)); label("Plasma (central) volume (mL/kg)") # workbook V_cen = 37.14

    # Minimal PBPK physiology (workbook)
    lvleaky <- fixed(log(62.4)); label("Leaky-tissue interstitial volume (mL/kg)") # workbook V_lea = 62.4
    lvtight <- fixed(log(115.9)); label("Tight-tissue interstitial volume (mL/kg)") # workbook V_tig = 115.9
    lvlymph <- fixed(log(44.4)); label("Lymph volume (mL/kg)") # workbook V_lym = 44.4
    llymphflow <- fixed(log(41.4288)); label("Total lymph flow (mL/day/kg)") # workbook L_D1 = 41.4288
    flymph_leaky <- fixed(0.65); label("Fraction of lymph flow through leaky tissue (fraction)") # workbook f_lea_D1 = 0.65
    flymph_tight <- fixed(0.35); label("Fraction of lymph flow through tight tissue (fraction)") # workbook f_tig_D1 = 0.35
    sigma_leaky <- fixed(0.67); label("Vascular reflection coefficient, leaky tissue (unitless)") # workbook sig_lea_D1 = 0.67
    sigma_tight <- fixed(0.95); label("Vascular reflection coefficient, tight tissue (unitless)") # workbook sig_tig_D1 = 0.95
    sigma_l <- fixed(0.2); label("Lymphatic reflection coefficient (unitless)") # workbook sig_lym_D1 = 0.2
    mw <- fixed(150000); label("Antibody molecular weight (g/mol)") # workbook MWab_D1 = 150000

    # Binding (R1 = CD20, R2 = CD3 in the published model)
    kon_cd20 <- fixed(3.6); label("Association rate constant for CD20 (1/nM/day)") # workbook kon_R1_D1 = 3.6
    kon_cd3 <- fixed(3.6); label("Association rate constant for CD3 (1/nM/day)") # workbook kon_R2_D1 = 3.6
    kd_cd20 <- fixed(68); label("Equilibrium dissociation constant for CD20 (nM)") # Table S2 KD,CD20 = 68 nM; workbook Kd_R1_D1 = 68
    kd_cd3 <- fixed(40); label("Equilibrium dissociation constant for CD3 (nM)") # Table S2 KD,CD3 = 40 nM; workbook Kd_R2_D1 = 40

    # Receptor capacity and turnover (central compartment only)
    bl_cd20 <- 4.4; label("Total CD20 capacity in the central compartment (nM)") # Table S2 R_CD20,central = 4.4 nM (fit; +/-30% range 3.3-5.8)
    bl_cd3 <- 19; label("Total CD3 capacity in the central compartment (nM)") # Table S2 R_CD3,central = 19 nM (fit; +/-30% range 13-25)
    thalf_cd20 <- fixed(1.0); label("CD20 turnover half-life (day)") # Table S2 T1/2,CD20 = 1.0 day; workbook Thalf_mR1 = 1
    thalf_cd3 <- fixed(0.7083); label("CD3 turnover half-life (day)") # workbook Thalf_mR2 = 0.7083 day; Table S2 T1/2,CD3 = 0.7 day

    propSd <- fixed(0); label("Proportional residual error (fraction; not estimated in the source)")
  })

  model({
    cl <- exp(lcl)
    vc <- exp(lvc)
    vleaky <- exp(lvleaky)
    vtight <- exp(lvtight)
    vlymph <- exp(lvlymph)
    lymphflow <- exp(llymphflow)
    l_leaky <- flymph_leaky * lymphflow
    l_tight <- flymph_tight * lymphflow

    k_c_leaky <- l_leaky * (1 - sigma_leaky) / vc
    k_c_tight <- l_tight * (1 - sigma_tight) / vc
    k_leaky_l <- l_leaky * (1 - sigma_l) / vleaky
    k_tight_l <- l_tight * (1 - sigma_l) / vtight
    k_l_c <- lymphflow / vlymph
    kel <- cl / vc

    # Receptor-bound drug internalizes at the receptor turnover rate (the
    # workbook's fold change in internalization, mR1_fcint / mR2_fcint, is 1).
    kdeg_cd20 <- log(2) / thalf_cd20
    kdeg_cd3 <- log(2) / thalf_cd3
    ksyn_cd20 <- bl_cd20 * kdeg_cd20
    ksyn_cd3 <- bl_cd3 * kdeg_cd3

    target_cd20_central(0) <- bl_cd20
    target_cd3_central(0) <- bl_cd3

    # Binding fluxes (nM/day). The second (trans) binding step has an
    # effective avidity of 1 (Data S1: trans-binding avidity constant = 1).
    b20 <- kon_cd20 * (central * target_cd20_central - kd_cd20 * drug_cd20_central)
    b3 <- kon_cd3 * (central * target_cd3_central - kd_cd3 * drug_cd3_central)
    b20to3 <- kon_cd3 * (drug_cd20_central * target_cd3_central - kd_cd3 * trimer)
    b3to20 <- kon_cd20 * (drug_cd3_central * target_cd20_central - kd_cd20 * trimer)

    d/dt(central) <- -kel * central - (k_c_leaky + k_c_tight) * central +
      k_l_c * vlymph / vc * lymph - b20 - b3
    d/dt(leaky) <- k_c_leaky * vc / vleaky * central - k_leaky_l * leaky
    d/dt(tight) <- k_c_tight * vc / vtight * central - k_tight_l * tight
    d/dt(lymph) <- k_leaky_l * vleaky / vlymph * leaky + k_tight_l * vtight / vlymph * tight - k_l_c * lymph

    d/dt(target_cd20_central) <- ksyn_cd20 - kdeg_cd20 * target_cd20_central - b20 - b3to20 +
      kdeg_cd3 * trimer
    d/dt(target_cd3_central) <- ksyn_cd3 - kdeg_cd3 * target_cd3_central - b3 - b20to3 +
      kdeg_cd20 * trimer
    d/dt(drug_cd20_central) <- b20 - b20to3 - kdeg_cd20 * drug_cd20_central
    d/dt(drug_cd3_central) <- b3 - b3to20 - kdeg_cd3 * drug_cd3_central
    d/dt(trimer) <- b20to3 + b3to20 - (kdeg_cd20 + kdeg_cd3) * trimer

    # Cumulative elimination by route (nmol/kg; Data S1 CL_D1_NS_nmol and
    # CL_D1_TMDD_R1/R2_cen_nmol without the body-weight factor)
    d/dt(elim_ns) <- kel * central * vc / 1000
    d/dt(elim_cd20) <- kdeg_cd20 * (drug_cd20_central + trimer) * vc / 1000
    d/dt(elim_cd3) <- kdeg_cd3 * (drug_cd3_central + trimer) * vc / 1000

    f(central) <- 1e9 / (mw * vc)

    # Free (unbound) drug in serum, the analyte of the cynomolgus PK assay
    Cc <- central * mw / 1e6
    total_drug_central <- central + drug_cd20_central + drug_cd3_central + trimer
    ro_cd20_central <- 100 * (1 - target_cd20_central / bl_cd20)
    ro_cd3_central <- 100 * (1 - target_cd3_central / bl_cd3)

    Cc ~ prop(propSd)
  })
}
