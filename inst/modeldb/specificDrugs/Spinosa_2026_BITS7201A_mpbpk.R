Spinosa_2026_BITS7201A_mpbpk <- function() {
  description <- paste(
    "QSP / minimal-PBPK. Generalized bispecific-antibody minimal PBPK-PD model",
    "(Spinosa 2026, Case Study 1) parameterized for BITS7201A, a bispecific",
    "antibody that neutralizes the two soluble cytokines IL-13 and IL-17AA, in",
    "healthy volunteers. Second-generation Cao-Jusko minimal PBPK backbone",
    "(plasma, leaky and tight interstitial tissue, lymph) with first-order",
    "subcutaneous absorption and nonspecific clearance from plasma only. Both",
    "cytokines are synthesized (zero order) and degraded (first order) in the",
    "plasma, leaky and tight compartments at the same baseline; the drug binds",
    "each cytokine reversibly and the single-bound complexes bind the second",
    "cytokine to form a drug-IL-13-IL-17AA ternary complex. All drug-cytokine",
    "complexes are cleared from plasma at the free-drug clearance and traffic",
    "between compartments at the free-drug rates, which is what makes total",
    "serum IL-17AA and IL-13 rise many-fold after dosing (target stabilization)."
  )
  reference <- paste(
    "Spinosa P, Joslyn L, Ramanujan S, Gadkar K, Hosseini I. A Generalized",
    "Minimal PBPK-PD Model of Bispecific Antibodies: Case Studies and",
    "Applications in Drug Development. CPT Pharmacometrics Syst Pharmacol.",
    "2026;15(2):e70167. doi:10.1002/psp4.70167. Model equations and default",
    "parameters from the Supporting Information Data S1 (SimBiology project",
    "PSP-2025-0197-s04.sbproj, model 'm1'; parameter workbook",
    "PSP-2025-0197-s07.xlsx sheet 'Case Study 1 Anti-IL17Anti-IL13'; Table S1).",
    "Clinical PK/PD data from Staton TL et al. BMC Pulm Med. 2019;19:5",
    "doi:10.1186/s12890-018-0763-9. The minimal PBPK backbone follows Cao Y,",
    "Balthasar JP, Jusko WJ. J Pharmacokinet Pharmacodyn. 2013;40:597-607",
    "(see modellib('Cao_2013_MEDI528')).",
    sep = " "
  )
  vignette <- "Spinosa_2026_bispecific_mpbpk"

  # Every state except `depot` holds a CONCENTRATION (nM) in a body-weight
  # normalised volume (mL/kg), as in the published SimBiology model, because
  # every binding rate law is written in nM. `depot` holds the subcutaneous
  # dose in mg. f(central) and the depot-to-central flux convert mg into nM:
  # mg / (mL/kg * kg) is g/L; dividing by the molecular weight (g/mol) gives
  # mol/L and multiplying by 1e9 gives nM.
  dosing <- c("depot", "central")
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  # The cytokines and the drug-cytokine complexes are located per mPBPK
  # compartment; the drug_<target>_<location> shape follows the registered
  # drug_cd3_central family, and trimer_<location> is the located form of the
  # registered `trimer` (drug bound to both of its targets).
  paper_specific_compartments <- c(
    "il13_central",
    "il13_leaky",
    "il13_tight",
    "il17_central",
    "il17_leaky",
    "il17_tight",
    "drug_il13_central",
    "drug_il13_leaky",
    "drug_il13_tight",
    "drug_il13_lymph",
    "drug_il17_central",
    "drug_il17_leaky",
    "drug_il17_tight",
    "drug_il17_lymph",
    "trimer_central",
    "trimer_leaky",
    "trimer_tight",
    "trimer_lymph"
  )

  compartmentData <- list(
    depot = list(analyte = "BITS7201A", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "BITS7201A (free)", units = "nM", specimen = "serum", verified = TRUE),
    leaky = list(analyte = "BITS7201A (free)", units = "nM", specimen = "tissue", verified = TRUE),
    tight = list(analyte = "BITS7201A (free)", units = "nM", specimen = "tissue", verified = TRUE),
    lymph = list(analyte = "BITS7201A (free)", units = "nM", specimen = "lymph", verified = TRUE),
    il13_central = list(analyte = "IL-13 (free)", units = "nM", specimen = "serum", verified = TRUE),
    il13_leaky = list(analyte = "IL-13 (free)", units = "nM", specimen = "tissue", verified = TRUE),
    il13_tight = list(analyte = "IL-13 (free)", units = "nM", specimen = "tissue", verified = TRUE),
    il17_central = list(analyte = "IL-17AA (free)", units = "nM", specimen = "serum", verified = TRUE),
    il17_leaky = list(analyte = "IL-17AA (free)", units = "nM", specimen = "tissue", verified = TRUE),
    il17_tight = list(analyte = "IL-17AA (free)", units = "nM", specimen = "tissue", verified = TRUE),
    drug_il13_central = list(analyte = "BITS7201A-IL-13 complex", units = "nM", specimen = "serum", verified = TRUE),
    drug_il13_leaky = list(analyte = "BITS7201A-IL-13 complex", units = "nM", specimen = "tissue", verified = TRUE),
    drug_il13_tight = list(analyte = "BITS7201A-IL-13 complex", units = "nM", specimen = "tissue", verified = TRUE),
    drug_il13_lymph = list(analyte = "BITS7201A-IL-13 complex", units = "nM", specimen = "lymph", verified = TRUE),
    drug_il17_central = list(analyte = "BITS7201A-IL-17AA complex", units = "nM", specimen = "serum", verified = TRUE),
    drug_il17_leaky = list(analyte = "BITS7201A-IL-17AA complex", units = "nM", specimen = "tissue", verified = TRUE),
    drug_il17_tight = list(analyte = "BITS7201A-IL-17AA complex", units = "nM", specimen = "tissue", verified = TRUE),
    drug_il17_lymph = list(analyte = "BITS7201A-IL-17AA complex", units = "nM", specimen = "lymph", verified = TRUE),
    trimer_central = list(
      analyte = "BITS7201A-IL-13-IL-17AA complex",
      units = "nM",
      specimen = "serum",
      verified = TRUE
    ),
    trimer_leaky = list(
      analyte = "BITS7201A-IL-13-IL-17AA complex",
      units = "nM",
      specimen = "tissue",
      verified = TRUE
    ),
    trimer_tight = list(
      analyte = "BITS7201A-IL-13-IL-17AA complex",
      units = "nM",
      specimen = "tissue",
      verified = TRUE
    ),
    trimer_lymph = list(analyte = "BITS7201A-IL-13-IL-17AA complex", units = "nM", specimen = "lymph", verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species = "human",
    n_subjects = NA_integer_,
    n_studies = 1L,
    age_range = "adults (healthy volunteers); not tabulated in Spinosa 2026",
    weight_range = "70 kg reference body weight (parameter workbook, Case Study 1)",
    sex_female_pct = NA_real_,
    race_ethnicity = "not reported",
    disease_state = "healthy volunteers (phase I single- and multiple-ascending-dose study of BITS7201A; Staton 2019)",
    dose_range = paste(
      "SAD: 30, 90 and 300 mg SC and 300 and 750 mg IV (model calibration);",
      "MAD: 150, 300 and 600 mg SC every 4 weeks x 3 (model validation)"
    ),
    regions = "not reported",
    notes = paste(
      "A mechanistic model calibrated to mean-level clinical PK and serum IL-17AA data;",
      "no subject-level fit, so there is no inter-individual variability and no",
      "residual-error estimate. Anti-drug-antibody positive and negative subjects",
      "were pooled (no apparent exposure effect). Baseline serum IL-17AA is low",
      "(often undetectable) and variable; the paper brackets it by simulating",
      "baselines across 0.04-0.3 pg/mL."
    )
  )

  ini({
    # Drug disposition (Table S1; workbook sheet 'Case Study 1 Anti-IL17Anti-IL13')
    lcl <- log(3.1); label("Nonspecific clearance from plasma (mL/day/kg)") # Table S1 CLns = 3.1 mL/kg/day (calibrated); workbook CLp_D1 = 3.1
    lvc <- log(52); label("Plasma (central) volume (mL/kg)") # Table S1 Vcen = 52 mL/kg (inferred from reported Cmax); workbook V_cen = 52
    lka <- log(0.2); label("Subcutaneous absorption rate constant (1/day)") # Table S1 kabs = 0.2 1/day; workbook kabs_D1 = 0.2
    lfdepot <- log(0.92); label("Subcutaneous bioavailability (fraction)") # Table S1 F = 0.92; workbook F_D1 = 0.92

    # Minimal PBPK physiology (workbook; Cao 2013 second-generation mPBPK values)
    lvleaky <- fixed(log(62.4)); label("Leaky-tissue interstitial volume (mL/kg)") # workbook V_lea = 62.4
    lvtight <- fixed(log(115.9)); label("Tight-tissue interstitial volume (mL/kg)") # workbook V_tig = 115.9
    lvlymph <- fixed(log(37.14)); label("Lymph volume (mL/kg)") # workbook V_lym = 37.14
    llymphflow <- fixed(log(41.4288)); label("Total lymph flow (mL/day/kg)") # workbook L_D1 = 41.4288
    flymph_leaky <- fixed(0.65); label("Fraction of lymph flow through leaky tissue (fraction)") # workbook f_lea_D1 = 0.65
    flymph_tight <- fixed(0.35); label("Fraction of lymph flow through tight tissue (fraction)") # workbook f_tig_D1 = 0.35
    sigma_leaky <- fixed(0.7); label("Vascular reflection coefficient, leaky tissue (unitless)") # workbook sig_lea_D1 = 0.7
    sigma_tight <- fixed(0.95); label("Vascular reflection coefficient, tight tissue (unitless)") # workbook sig_tig_D1 = 0.95
    sigma_l <- fixed(0.2); label("Lymphatic reflection coefficient (unitless)") # workbook sig_lym_D1 = 0.2
    mw <- fixed(150000); label("Antibody molecular weight (g/mol)") # workbook MWab_D1 = 150000; Data S1 text '150 kDa'
    bw <- fixed(70); label("Body weight used to scale the per-kg volumes for flat-mg doses (kg)") # workbook BW = 70

    # Binding (R1 = IL-13, R2 = IL-17AA in the published model)
    kon_il13 <- fixed(3.6); label("Association rate constant for IL-13 (1/nM/day)") # workbook kon_R1_D1 = 3.6
    kon_il17 <- fixed(3.6); label("Association rate constant for IL-17AA (1/nM/day)") # workbook kon_R2_D1 = 3.6
    kd_il13 <- fixed(0.005); label("Equilibrium dissociation constant for IL-13 (nM)") # Table S1 KD = 0.005 nM (IL-13); workbook Kd_R1_D1 = 0.005
    kd_il17 <- fixed(0.010); label("Equilibrium dissociation constant for IL-17AA (nM)") # Table S1 KD = 0.010 nM (IL-17AA); workbook Kd_R2_D1 = 0.01

    # Cytokine turnover (baseline identical in plasma, leaky and tight tissue)
    bl_il13 <- 6e-6; label("Baseline free IL-13 in plasma and tissue (nM)") # workbook sR1_0_cen/lea/tig = 6e-06 nM (Table S1 Css 0.01-0.1 pg/mL)
    bl_il17 <- 2.83e-6; label("Baseline free IL-17AA in plasma and tissue (nM)") # workbook sR2_0_cen/lea/tig = 2.83e-06 nM (fitted; Table S1 Css range)
    thalf_il13 <- fixed(0.0208); label("IL-13 turnover half-life (day)") # workbook Thalf_sR1 = 0.0208 day; Table S1 T1/2 = 30 min
    thalf_il17 <- 0.0056; label("IL-17AA turnover half-life (day)") # workbook Thalf_sR2 = 0.0056 day (fitted); Table S1 T1/2 = 8 min

    propSd <- fixed(0); label("Proportional residual error (fraction; not estimated in the source)")
  })

  model({
    cl <- exp(lcl)
    vc <- exp(lvc)
    ka <- exp(lka)
    vleaky <- exp(lvleaky)
    vtight <- exp(lvtight)
    vlymph <- exp(lvlymph)
    lymphflow <- exp(llymphflow)
    l_leaky <- flymph_leaky * lymphflow
    l_tight <- flymph_tight * lymphflow

    # Plasma-to-tissue, tissue-to-lymph and lymph-to-plasma first-order rate
    # constants (1/day) in the units of the compartment being left; the
    # receiving compartment gains the same amount, i.e. the flux times the
    # volume ratio (Data S1 'transport_between_compartments').
    k_c_leaky <- l_leaky * (1 - sigma_leaky) / vc
    k_c_tight <- l_tight * (1 - sigma_tight) / vc
    k_leaky_l <- l_leaky * (1 - sigma_l) / vleaky
    k_tight_l <- l_tight * (1 - sigma_l) / vtight
    k_l_c <- lymphflow / vlymph
    kel <- cl / vc

    kdeg_il13 <- log(2) / thalf_il13
    kdeg_il17 <- log(2) / thalf_il17
    ksyn_il13 <- bl_il13 * kdeg_il13
    ksyn_il17 <- bl_il17 * kdeg_il17

    # mg -> nM in plasma
    conv <- 1e9 / (mw * vc * bw)

    il13_central(0) <- bl_il13
    il13_leaky(0) <- bl_il13
    il13_tight(0) <- bl_il13
    il17_central(0) <- bl_il17
    il17_leaky(0) <- bl_il17
    il17_tight(0) <- bl_il17

    # Binding fluxes (nM/day): drug + cytokine, and single-bound complex +
    # the second cytokine forming the ternary complex.
    b13_c <- kon_il13 * (central * il13_central - kd_il13 * drug_il13_central)
    b17_c <- kon_il17 * (central * il17_central - kd_il17 * drug_il17_central)
    b17to13_c <- kon_il13 * (drug_il17_central * il13_central - kd_il13 * trimer_central)
    b13to17_c <- kon_il17 * (drug_il13_central * il17_central - kd_il17 * trimer_central)

    b13_lk <- kon_il13 * (leaky * il13_leaky - kd_il13 * drug_il13_leaky)
    b17_lk <- kon_il17 * (leaky * il17_leaky - kd_il17 * drug_il17_leaky)
    b17to13_lk <- kon_il13 * (drug_il17_leaky * il13_leaky - kd_il13 * trimer_leaky)
    b13to17_lk <- kon_il17 * (drug_il13_leaky * il17_leaky - kd_il17 * trimer_leaky)

    b13_tt <- kon_il13 * (tight * il13_tight - kd_il13 * drug_il13_tight)
    b17_tt <- kon_il17 * (tight * il17_tight - kd_il17 * drug_il17_tight)
    b17to13_tt <- kon_il13 * (drug_il17_tight * il13_tight - kd_il13 * trimer_tight)
    b13to17_tt <- kon_il17 * (drug_il13_tight * il17_tight - kd_il17 * trimer_tight)

    # Free drug
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot * conv - kel * central - (k_c_leaky + k_c_tight) * central +
      k_l_c * vlymph / vc * lymph - b13_c - b17_c
    d/dt(leaky) <- k_c_leaky * vc / vleaky * central - k_leaky_l * leaky - b13_lk - b17_lk
    d/dt(tight) <- k_c_tight * vc / vtight * central - k_tight_l * tight - b13_tt - b17_tt
    d/dt(lymph) <- k_leaky_l * vleaky / vlymph * leaky + k_tight_l * vtight / vlymph * tight - k_l_c * lymph

    # Free cytokines (no trafficking of free cytokine)
    d/dt(il13_central) <- ksyn_il13 - kdeg_il13 * il13_central - b13_c - b17to13_c
    d/dt(il13_leaky) <- ksyn_il13 - kdeg_il13 * il13_leaky - b13_lk - b17to13_lk
    d/dt(il13_tight) <- ksyn_il13 - kdeg_il13 * il13_tight - b13_tt - b17to13_tt
    d/dt(il17_central) <- ksyn_il17 - kdeg_il17 * il17_central - b17_c - b13to17_c
    d/dt(il17_leaky) <- ksyn_il17 - kdeg_il17 * il17_leaky - b17_lk - b13to17_lk
    d/dt(il17_tight) <- ksyn_il17 - kdeg_il17 * il17_tight - b17_tt - b13to17_tt

    # Drug-cytokine complexes: cleared from plasma at the free-drug clearance
    # and trafficking at the free-drug rates
    d/dt(drug_il13_central) <- b13_c - b13to17_c - kel * drug_il13_central -
      (k_c_leaky + k_c_tight) * drug_il13_central + k_l_c * vlymph / vc * drug_il13_lymph
    d/dt(drug_il13_leaky) <- b13_lk - b13to17_lk + k_c_leaky * vc / vleaky * drug_il13_central -
      k_leaky_l * drug_il13_leaky
    d/dt(drug_il13_tight) <- b13_tt - b13to17_tt + k_c_tight * vc / vtight * drug_il13_central -
      k_tight_l * drug_il13_tight
    d/dt(drug_il13_lymph) <- k_leaky_l * vleaky / vlymph * drug_il13_leaky +
      k_tight_l * vtight / vlymph * drug_il13_tight - k_l_c * drug_il13_lymph

    d/dt(drug_il17_central) <- b17_c - b17to13_c - kel * drug_il17_central -
      (k_c_leaky + k_c_tight) * drug_il17_central + k_l_c * vlymph / vc * drug_il17_lymph
    d/dt(drug_il17_leaky) <- b17_lk - b17to13_lk + k_c_leaky * vc / vleaky * drug_il17_central -
      k_leaky_l * drug_il17_leaky
    d/dt(drug_il17_tight) <- b17_tt - b17to13_tt + k_c_tight * vc / vtight * drug_il17_central -
      k_tight_l * drug_il17_tight
    d/dt(drug_il17_lymph) <- k_leaky_l * vleaky / vlymph * drug_il17_leaky +
      k_tight_l * vtight / vlymph * drug_il17_tight - k_l_c * drug_il17_lymph

    d/dt(trimer_central) <- b17to13_c + b13to17_c - kel * trimer_central -
      (k_c_leaky + k_c_tight) * trimer_central + k_l_c * vlymph / vc * trimer_lymph
    d/dt(trimer_leaky) <- b17to13_lk + b13to17_lk + k_c_leaky * vc / vleaky * trimer_central -
      k_leaky_l * trimer_leaky
    d/dt(trimer_tight) <- b17to13_tt + b13to17_tt + k_c_tight * vc / vtight * trimer_central -
      k_tight_l * trimer_tight
    d/dt(trimer_lymph) <- k_leaky_l * vleaky / vlymph * trimer_leaky +
      k_tight_l * vtight / vlymph * trimer_tight - k_l_c * trimer_lymph

    f(depot) <- exp(lfdepot)
    f(central) <- conv

    # Total (free + cytokine-bound) drug in serum, as measured by the PK assay
    total_drug_central <- central + drug_il13_central + drug_il17_central + trimer_central
    Cc <- total_drug_central * mw / 1e6

    # Total serum cytokines (free + bound; nM) and their fold over baseline
    il17_total_central <- il17_central + drug_il17_central + trimer_central
    il13_total_central <- il13_central + drug_il13_central + trimer_central
    il17_fold_central <- il17_total_central / bl_il17
    il13_fold_central <- il13_total_central / bl_il13

    # Free cytokine as percent of baseline (Figure 3d,e 'target neutralization')
    il17_pctfree_central <- 100 * il17_central / bl_il17
    il17_pctfree_leaky <- 100 * il17_leaky / bl_il17
    il17_pctfree_tight <- 100 * il17_tight / bl_il17
    il13_pctfree_central <- 100 * il13_central / bl_il13
    il13_pctfree_leaky <- 100 * il13_leaky / bl_il13
    il13_pctfree_tight <- 100 * il13_tight / bl_il13

    Cc ~ prop(propSd)
  })
}
