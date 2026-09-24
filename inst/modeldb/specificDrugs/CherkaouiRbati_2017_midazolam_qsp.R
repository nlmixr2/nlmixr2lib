CherkaouiRbati_2017_midazolam_qsp <- function() {
  description <- paste(
    "QSP / PBPK (7 compartments plus a spatially resolved liver lobule,",
    "122 ODE states). Dynamic CYP3A4 drug-drug-interaction model in which",
    "the liver is not well stirred: the Fig 1C lobule algorithm generates",
    "5 sinusoid levels whose radius narrows and blood velocity rises",
    "toward the central vein, and the resulting 1-D convection equation is",
    "solved by the method of lines over 20 intervals (21 nodes,",
    "`sinusoid_slab<n>` / `hepatocyte_slab<n>`). Drug crosses the",
    "sinusoidal membrane by passive permeability and is metabolised inside",
    "the hepatocytes, so enzyme level, inhibitor concentration and",
    "metabolic rate all vary with position along the sinusoid (paper Fig",
    "7). Two drugs are carried simultaneously - midazolam as the CYP3A4",
    "victim probe and one perpetrator - and all three interaction",
    "mechanisms act on a shared, spatially resolved CYP3A4 pool:",
    "competitive inhibition, mechanism-based inactivation and additive",
    "induction. A two-sub-compartment gut (enterocytes plus portal vein)",
    "supplies first-pass metabolism. The perpetrator slot is",
    "parameterised, and ships set to ketoconazole; the validation vignette",
    "swaps in each of the paper's other perpetrators to reproduce the",
    "Table 9 interaction ratios. Deterministic: no IIV and no residual",
    "error are reported or implemented.",
    sep = " "
  )

  reference <- paste(
    "Cherkaoui-Rbati MH, Paine SW, Littlewood P, Rauch C. A quantitative",
    "systems pharmacology approach, incorporating a novel liver model, for",
    "predicting pharmacokinetic drug-drug interactions. PLoS One.",
    "2017;12(9):e0183794. doi:10.1371/journal.pone.0183794",
    sep = " "
  )

  vignette <- "CherkaouiRbati_2017_midazolam_qsp"

  units <- list(time = "h", dosing = "umol", concentration = "umol/L")

  # This model carries no covariates: every parameter is a fixed physiological
  # or in vitro constant for an average 70 kg man (S1 Table), and the paper
  # reports no covariate model.
  covariateData <- list()

  compartmentData <- list(
    depot = list(analyte = "midazolam", units = "umol", specimen = "administration site", verified = TRUE),
    depot_perpetrator = list(
      analyte = "ketoconazole",
      units = "umol",
      specimen = "administration site",
      verified = TRUE
    ),
    a_gut = list(analyte = "midazolam", units = "umol", specimen = "tissue", verified = TRUE),
    a_gut_perpetrator = list(analyte = "ketoconazole", units = "umol", specimen = "tissue", verified = TRUE),
    a_portal = list(analyte = "midazolam", units = "umol", specimen = "whole blood", verified = TRUE),
    a_portal_perpetrator = list(analyte = "ketoconazole", units = "umol", specimen = "whole blood", verified = TRUE),
    a_arterial = list(analyte = "midazolam", units = "umol", specimen = "whole blood", verified = TRUE),
    a_arterial_perpetrator = list(analyte = "ketoconazole", units = "umol", specimen = "whole blood", verified = TRUE),
    a_venous = list(analyte = "midazolam", units = "umol", specimen = "whole blood", verified = TRUE),
    a_venous_perpetrator = list(analyte = "ketoconazole", units = "umol", specimen = "whole blood", verified = TRUE),
    a_kidney = list(analyte = "midazolam", units = "umol", specimen = "tissue", verified = TRUE),
    a_kidney_perpetrator = list(analyte = "ketoconazole", units = "umol", specimen = "tissue", verified = TRUE),
    a_lung = list(analyte = "midazolam", units = "umol", specimen = "tissue", verified = TRUE),
    a_lung_perpetrator = list(analyte = "ketoconazole", units = "umol", specimen = "tissue", verified = TRUE),
    a_remainder = list(analyte = "midazolam", units = "umol", specimen = "tissue", verified = TRUE),
    a_remainder_perpetrator = list(analyte = "ketoconazole", units = "umol", specimen = "tissue", verified = TRUE),
    enzyme_gut = list(
      analyte = "cytochrome P450 3A4",
      units = "fold of baseline",
      specimen = "tissue",
      verified = TRUE
    ),
    sinusoid_slab1 = list(analyte = "midazolam", units = "umol", specimen = "whole blood", verified = TRUE),
    sinusoid_perpetrator_slab1 = list(
      analyte = "ketoconazole",
      units = "umol",
      specimen = "whole blood",
      verified = TRUE
    ),
    hepatocyte_slab1 = list(analyte = "midazolam", units = "umol", specimen = "tissue", verified = TRUE),
    hepatocyte_perpetrator_slab1 = list(analyte = "ketoconazole", units = "umol", specimen = "tissue", verified = TRUE),
    enzyme_liver_slab1 = list(
      analyte = "cytochrome P450 3A4",
      units = "fold of baseline",
      specimen = "tissue",
      verified = TRUE
    ),
    sinusoid_slab2 = list(analyte = "midazolam", units = "umol", specimen = "whole blood", verified = TRUE),
    sinusoid_perpetrator_slab2 = list(
      analyte = "ketoconazole",
      units = "umol",
      specimen = "whole blood",
      verified = TRUE
    ),
    hepatocyte_slab2 = list(analyte = "midazolam", units = "umol", specimen = "tissue", verified = TRUE),
    hepatocyte_perpetrator_slab2 = list(analyte = "ketoconazole", units = "umol", specimen = "tissue", verified = TRUE),
    enzyme_liver_slab2 = list(
      analyte = "cytochrome P450 3A4",
      units = "fold of baseline",
      specimen = "tissue",
      verified = TRUE
    ),
    sinusoid_slab3 = list(analyte = "midazolam", units = "umol", specimen = "whole blood", verified = TRUE),
    sinusoid_perpetrator_slab3 = list(
      analyte = "ketoconazole",
      units = "umol",
      specimen = "whole blood",
      verified = TRUE
    ),
    hepatocyte_slab3 = list(analyte = "midazolam", units = "umol", specimen = "tissue", verified = TRUE),
    hepatocyte_perpetrator_slab3 = list(analyte = "ketoconazole", units = "umol", specimen = "tissue", verified = TRUE),
    enzyme_liver_slab3 = list(
      analyte = "cytochrome P450 3A4",
      units = "fold of baseline",
      specimen = "tissue",
      verified = TRUE
    ),
    sinusoid_slab4 = list(analyte = "midazolam", units = "umol", specimen = "whole blood", verified = TRUE),
    sinusoid_perpetrator_slab4 = list(
      analyte = "ketoconazole",
      units = "umol",
      specimen = "whole blood",
      verified = TRUE
    ),
    hepatocyte_slab4 = list(analyte = "midazolam", units = "umol", specimen = "tissue", verified = TRUE),
    hepatocyte_perpetrator_slab4 = list(analyte = "ketoconazole", units = "umol", specimen = "tissue", verified = TRUE),
    enzyme_liver_slab4 = list(
      analyte = "cytochrome P450 3A4",
      units = "fold of baseline",
      specimen = "tissue",
      verified = TRUE
    ),
    sinusoid_slab5 = list(analyte = "midazolam", units = "umol", specimen = "whole blood", verified = TRUE),
    sinusoid_perpetrator_slab5 = list(
      analyte = "ketoconazole",
      units = "umol",
      specimen = "whole blood",
      verified = TRUE
    ),
    hepatocyte_slab5 = list(analyte = "midazolam", units = "umol", specimen = "tissue", verified = TRUE),
    hepatocyte_perpetrator_slab5 = list(analyte = "ketoconazole", units = "umol", specimen = "tissue", verified = TRUE),
    enzyme_liver_slab5 = list(
      analyte = "cytochrome P450 3A4",
      units = "fold of baseline",
      specimen = "tissue",
      verified = TRUE
    ),
    sinusoid_slab6 = list(analyte = "midazolam", units = "umol", specimen = "whole blood", verified = TRUE),
    sinusoid_perpetrator_slab6 = list(
      analyte = "ketoconazole",
      units = "umol",
      specimen = "whole blood",
      verified = TRUE
    ),
    hepatocyte_slab6 = list(analyte = "midazolam", units = "umol", specimen = "tissue", verified = TRUE),
    hepatocyte_perpetrator_slab6 = list(analyte = "ketoconazole", units = "umol", specimen = "tissue", verified = TRUE),
    enzyme_liver_slab6 = list(
      analyte = "cytochrome P450 3A4",
      units = "fold of baseline",
      specimen = "tissue",
      verified = TRUE
    ),
    sinusoid_slab7 = list(analyte = "midazolam", units = "umol", specimen = "whole blood", verified = TRUE),
    sinusoid_perpetrator_slab7 = list(
      analyte = "ketoconazole",
      units = "umol",
      specimen = "whole blood",
      verified = TRUE
    ),
    hepatocyte_slab7 = list(analyte = "midazolam", units = "umol", specimen = "tissue", verified = TRUE),
    hepatocyte_perpetrator_slab7 = list(analyte = "ketoconazole", units = "umol", specimen = "tissue", verified = TRUE),
    enzyme_liver_slab7 = list(
      analyte = "cytochrome P450 3A4",
      units = "fold of baseline",
      specimen = "tissue",
      verified = TRUE
    ),
    sinusoid_slab8 = list(analyte = "midazolam", units = "umol", specimen = "whole blood", verified = TRUE),
    sinusoid_perpetrator_slab8 = list(
      analyte = "ketoconazole",
      units = "umol",
      specimen = "whole blood",
      verified = TRUE
    ),
    hepatocyte_slab8 = list(analyte = "midazolam", units = "umol", specimen = "tissue", verified = TRUE),
    hepatocyte_perpetrator_slab8 = list(analyte = "ketoconazole", units = "umol", specimen = "tissue", verified = TRUE),
    enzyme_liver_slab8 = list(
      analyte = "cytochrome P450 3A4",
      units = "fold of baseline",
      specimen = "tissue",
      verified = TRUE
    ),
    sinusoid_slab9 = list(analyte = "midazolam", units = "umol", specimen = "whole blood", verified = TRUE),
    sinusoid_perpetrator_slab9 = list(
      analyte = "ketoconazole",
      units = "umol",
      specimen = "whole blood",
      verified = TRUE
    ),
    hepatocyte_slab9 = list(analyte = "midazolam", units = "umol", specimen = "tissue", verified = TRUE),
    hepatocyte_perpetrator_slab9 = list(analyte = "ketoconazole", units = "umol", specimen = "tissue", verified = TRUE),
    enzyme_liver_slab9 = list(
      analyte = "cytochrome P450 3A4",
      units = "fold of baseline",
      specimen = "tissue",
      verified = TRUE
    ),
    sinusoid_slab10 = list(analyte = "midazolam", units = "umol", specimen = "whole blood", verified = TRUE),
    sinusoid_perpetrator_slab10 = list(
      analyte = "ketoconazole",
      units = "umol",
      specimen = "whole blood",
      verified = TRUE
    ),
    hepatocyte_slab10 = list(analyte = "midazolam", units = "umol", specimen = "tissue", verified = TRUE),
    hepatocyte_perpetrator_slab10 = list(
      analyte = "ketoconazole",
      units = "umol",
      specimen = "tissue",
      verified = TRUE
    ),
    enzyme_liver_slab10 = list(
      analyte = "cytochrome P450 3A4",
      units = "fold of baseline",
      specimen = "tissue",
      verified = TRUE
    ),
    sinusoid_slab11 = list(analyte = "midazolam", units = "umol", specimen = "whole blood", verified = TRUE),
    sinusoid_perpetrator_slab11 = list(
      analyte = "ketoconazole",
      units = "umol",
      specimen = "whole blood",
      verified = TRUE
    ),
    hepatocyte_slab11 = list(analyte = "midazolam", units = "umol", specimen = "tissue", verified = TRUE),
    hepatocyte_perpetrator_slab11 = list(
      analyte = "ketoconazole",
      units = "umol",
      specimen = "tissue",
      verified = TRUE
    ),
    enzyme_liver_slab11 = list(
      analyte = "cytochrome P450 3A4",
      units = "fold of baseline",
      specimen = "tissue",
      verified = TRUE
    ),
    sinusoid_slab12 = list(analyte = "midazolam", units = "umol", specimen = "whole blood", verified = TRUE),
    sinusoid_perpetrator_slab12 = list(
      analyte = "ketoconazole",
      units = "umol",
      specimen = "whole blood",
      verified = TRUE
    ),
    hepatocyte_slab12 = list(analyte = "midazolam", units = "umol", specimen = "tissue", verified = TRUE),
    hepatocyte_perpetrator_slab12 = list(
      analyte = "ketoconazole",
      units = "umol",
      specimen = "tissue",
      verified = TRUE
    ),
    enzyme_liver_slab12 = list(
      analyte = "cytochrome P450 3A4",
      units = "fold of baseline",
      specimen = "tissue",
      verified = TRUE
    ),
    sinusoid_slab13 = list(analyte = "midazolam", units = "umol", specimen = "whole blood", verified = TRUE),
    sinusoid_perpetrator_slab13 = list(
      analyte = "ketoconazole",
      units = "umol",
      specimen = "whole blood",
      verified = TRUE
    ),
    hepatocyte_slab13 = list(analyte = "midazolam", units = "umol", specimen = "tissue", verified = TRUE),
    hepatocyte_perpetrator_slab13 = list(
      analyte = "ketoconazole",
      units = "umol",
      specimen = "tissue",
      verified = TRUE
    ),
    enzyme_liver_slab13 = list(
      analyte = "cytochrome P450 3A4",
      units = "fold of baseline",
      specimen = "tissue",
      verified = TRUE
    ),
    sinusoid_slab14 = list(analyte = "midazolam", units = "umol", specimen = "whole blood", verified = TRUE),
    sinusoid_perpetrator_slab14 = list(
      analyte = "ketoconazole",
      units = "umol",
      specimen = "whole blood",
      verified = TRUE
    ),
    hepatocyte_slab14 = list(analyte = "midazolam", units = "umol", specimen = "tissue", verified = TRUE),
    hepatocyte_perpetrator_slab14 = list(
      analyte = "ketoconazole",
      units = "umol",
      specimen = "tissue",
      verified = TRUE
    ),
    enzyme_liver_slab14 = list(
      analyte = "cytochrome P450 3A4",
      units = "fold of baseline",
      specimen = "tissue",
      verified = TRUE
    ),
    sinusoid_slab15 = list(analyte = "midazolam", units = "umol", specimen = "whole blood", verified = TRUE),
    sinusoid_perpetrator_slab15 = list(
      analyte = "ketoconazole",
      units = "umol",
      specimen = "whole blood",
      verified = TRUE
    ),
    hepatocyte_slab15 = list(analyte = "midazolam", units = "umol", specimen = "tissue", verified = TRUE),
    hepatocyte_perpetrator_slab15 = list(
      analyte = "ketoconazole",
      units = "umol",
      specimen = "tissue",
      verified = TRUE
    ),
    enzyme_liver_slab15 = list(
      analyte = "cytochrome P450 3A4",
      units = "fold of baseline",
      specimen = "tissue",
      verified = TRUE
    ),
    sinusoid_slab16 = list(analyte = "midazolam", units = "umol", specimen = "whole blood", verified = TRUE),
    sinusoid_perpetrator_slab16 = list(
      analyte = "ketoconazole",
      units = "umol",
      specimen = "whole blood",
      verified = TRUE
    ),
    hepatocyte_slab16 = list(analyte = "midazolam", units = "umol", specimen = "tissue", verified = TRUE),
    hepatocyte_perpetrator_slab16 = list(
      analyte = "ketoconazole",
      units = "umol",
      specimen = "tissue",
      verified = TRUE
    ),
    enzyme_liver_slab16 = list(
      analyte = "cytochrome P450 3A4",
      units = "fold of baseline",
      specimen = "tissue",
      verified = TRUE
    ),
    sinusoid_slab17 = list(analyte = "midazolam", units = "umol", specimen = "whole blood", verified = TRUE),
    sinusoid_perpetrator_slab17 = list(
      analyte = "ketoconazole",
      units = "umol",
      specimen = "whole blood",
      verified = TRUE
    ),
    hepatocyte_slab17 = list(analyte = "midazolam", units = "umol", specimen = "tissue", verified = TRUE),
    hepatocyte_perpetrator_slab17 = list(
      analyte = "ketoconazole",
      units = "umol",
      specimen = "tissue",
      verified = TRUE
    ),
    enzyme_liver_slab17 = list(
      analyte = "cytochrome P450 3A4",
      units = "fold of baseline",
      specimen = "tissue",
      verified = TRUE
    ),
    sinusoid_slab18 = list(analyte = "midazolam", units = "umol", specimen = "whole blood", verified = TRUE),
    sinusoid_perpetrator_slab18 = list(
      analyte = "ketoconazole",
      units = "umol",
      specimen = "whole blood",
      verified = TRUE
    ),
    hepatocyte_slab18 = list(analyte = "midazolam", units = "umol", specimen = "tissue", verified = TRUE),
    hepatocyte_perpetrator_slab18 = list(
      analyte = "ketoconazole",
      units = "umol",
      specimen = "tissue",
      verified = TRUE
    ),
    enzyme_liver_slab18 = list(
      analyte = "cytochrome P450 3A4",
      units = "fold of baseline",
      specimen = "tissue",
      verified = TRUE
    ),
    sinusoid_slab19 = list(analyte = "midazolam", units = "umol", specimen = "whole blood", verified = TRUE),
    sinusoid_perpetrator_slab19 = list(
      analyte = "ketoconazole",
      units = "umol",
      specimen = "whole blood",
      verified = TRUE
    ),
    hepatocyte_slab19 = list(analyte = "midazolam", units = "umol", specimen = "tissue", verified = TRUE),
    hepatocyte_perpetrator_slab19 = list(
      analyte = "ketoconazole",
      units = "umol",
      specimen = "tissue",
      verified = TRUE
    ),
    enzyme_liver_slab19 = list(
      analyte = "cytochrome P450 3A4",
      units = "fold of baseline",
      specimen = "tissue",
      verified = TRUE
    ),
    sinusoid_slab20 = list(analyte = "midazolam", units = "umol", specimen = "whole blood", verified = TRUE),
    sinusoid_perpetrator_slab20 = list(
      analyte = "ketoconazole",
      units = "umol",
      specimen = "whole blood",
      verified = TRUE
    ),
    hepatocyte_slab20 = list(analyte = "midazolam", units = "umol", specimen = "tissue", verified = TRUE),
    hepatocyte_perpetrator_slab20 = list(
      analyte = "ketoconazole",
      units = "umol",
      specimen = "tissue",
      verified = TRUE
    ),
    enzyme_liver_slab20 = list(
      analyte = "cytochrome P450 3A4",
      units = "fold of baseline",
      specimen = "tissue",
      verified = TRUE
    ),
    sinusoid_slab21 = list(analyte = "midazolam", units = "umol", specimen = "whole blood", verified = TRUE),
    sinusoid_perpetrator_slab21 = list(
      analyte = "ketoconazole",
      units = "umol",
      specimen = "whole blood",
      verified = TRUE
    ),
    hepatocyte_slab21 = list(analyte = "midazolam", units = "umol", specimen = "tissue", verified = TRUE),
    hepatocyte_perpetrator_slab21 = list(
      analyte = "ketoconazole",
      units = "umol",
      specimen = "tissue",
      verified = TRUE
    ),
    enzyme_liver_slab21 = list(
      analyte = "cytochrome P450 3A4",
      units = "fold of baseline",
      specimen = "tissue",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = NA_integer_,
    n_studies = 10L,
    disease_state = "healthy adult volunteers in the 10 published CYP3A4 interaction studies of Table 2",
    dose_range = "midazolam 5-15 mg orally; perpetrators 0.03-600 mg orally, single or repeated (Table 2)",
    notes = paste(
      "Not fitted to individual data. Every physiological volume and blood flow",
      "is the average value for a 70 kg man (S1 Table) and every drug parameter",
      "is taken from the literature or computed with published algorithms; the",
      "Discussion states that no parameter was fitted. The model was evaluated",
      "against the AUC and Cmax interaction ratios of 10 clinical studies",
      "(Table 2 and Table 9), all using midazolam as the CYP3A4 probe.",
      sep = " "
    )
  )

  ini({
    # ---- Physiology: S1 Table, average 70 kg man (flows L/min -> L/h) ----
    vArt <- fixed(1.73); label("Arterial blood volume (L)")  # S1 Table, arterial blood
    vVen <- fixed(3.47); label("Venous blood volume (L)")  # S1 Table, venous blood
    vKid <- fixed(0.28); label("Kidney volume (L)")  # S1 Table, kidney
    vLun <- fixed(1.17); label("Lung volume (L)")  # S1 Table, lungs
    vRem <- fixed(63.67); label("Rest-of-body volume (L)")  # S1 Table, sum of adipose, bone, brain, heart, muscle, pancreas, skin, spleen and stomach
    vGut <- fixed(0.23); label("Gut-wall volume (L)")  # S1 Table, gut wall
    vPor <- fixed(0.34); label("Portal vein volume (L)")  # S1 Table, portal vein
    qTot <- fixed(326.46); label("Cardiac output (L/h)")  # S1 Table, total blood flow 5.441 L/min
    qKid <- fixed(74.4); label("Kidney blood flow (L/h)")  # S1 Table, 1.24 L/min
    qLiv <- fixed(87); label("Total liver blood flow (L/h)")  # S1 Table, 1.45 L/min
    qPor <- fixed(66); label("Portal vein blood flow (L/h)")  # S1 Table, 1.1 L/min
    qHa <- fixed(21); label("Hepatic artery blood flow (L/h)")  # S1 Table, liver minus portal vein
    qRem <- fixed(165.06); label("Rest-of-body blood flow (L/h)")  # S1 Table, cardiac output minus liver and kidney

    # ---- CYP3A4: S2 Table amounts and degradation rates ----
    e0Liv <- fixed(6.62954); label("Baseline hepatic CYP3A4 concentration (umol/L)")  # S2 Table 9228.32 nmol / hepatocyte volume 1.392 L from the lobule model
    e0Gut <- fixed(0.304348); label("Baseline gut-wall CYP3A4 concentration (umol/L)")  # S2 Table 70.00 nmol / gut-wall volume 0.23 L
    kdegLiv <- fixed(0.0192); label("Hepatic CYP3A4 degradation rate (1/h)")  # S2 Table, k_deg
    kdegGut <- fixed(0.0288); label("Gut-wall CYP3A4 degradation rate (1/h)")  # S2 Table, k_deg^g

    # ---- Midazolam (victim) ----
    fubMdz <- fixed(0.04); label("Midazolam (victim) blood fraction unbound")  # Midazolam row of Table 4 column 'f_u^b'
    fuhMdz <- fixed(0.0202); label("Midazolam (victim) hepatocyte fraction unbound")  # Midazolam row of Table 4 column 'f_u^h'
    fugMdz <- fixed(0.0189); label("Midazolam (victim) gut-wall fraction unbound")  # Midazolam row of Table 4 column 'f_u^gw'
    rbpMdz <- fixed(0.66); label("Midazolam (victim) blood-to-plasma ratio")  # Midazolam row of Table 4 column 'R_BP'
    kpRemMdz <- fixed(0.84); label("Midazolam (victim) rest-of-body-to-plasma partition coefficient")  # Midazolam row of Table 5 column 'K_p,RB'
    kpKidMdz <- fixed(1.41); label("Midazolam (victim) kidney-to-plasma partition coefficient")  # Midazolam row of Table 5 column 'K_p,Kidney'
    kpLunMdz <- fixed(1.61); label("Midazolam (victim) lung-to-plasma partition coefficient")  # Midazolam row of Table 5 column 'K_p,Lungs'
    permMdz <- fixed(24228); label("Midazolam (victim) hepatocyte membrane permeability (um/h)")  # Midazolam row of Table 7 column 'P'
    km2Mdz <- fixed(2.3); label("Midazolam (victim) Michaelis constant of the non-CYP3A4 pathway K_m,2 (umol/L)")  # Midazolam row of Table 3 column 'K_m'
    rkm1Mdz <- fixed(0.434783); label("Midazolam (victim) reciprocal CYP3A4 Michaelis constant 1/K_m,1 (L/umol)")  # Midazolam row of Table 3 'K_m' with f_m,3A4 = 0.96; the drug is a CYP3A4 substrate so K_m,1 is the tabulated K_m
    cli3a4Mdz <- fixed(1373.24); label("Midazolam (victim) CYP3A4 intrinsic clearance per unit hepatocyte volume k_cat*E_0/K_m,1 (1/h)")  # Eq 21 line 1 as f_m,3A4 * CL*_m,int / V_h = 0.96 * 1991.2 / 1.392, Table 3 and Results V_h
    vmax2Mdz <- fixed(131.602); label("Midazolam (victim) hepatic non-CYP3A4 V_max,2 (umol/L/h)")  # Eq 21 line 2 as (1 - 0.96) * 1991.2 * 2.3 / 1.392
    vmax2gMdz <- fixed(0.0481081); label("Midazolam (victim) gut-wall non-CYP3A4 V_max,2^g (umol/L/h)")  # Eq 21 line 3 from S2 Table CYP fractions and amounts
    rkiMdz <- fixed(0); label("Midazolam (victim) reciprocal competitive inhibition constant 1/K_i (L/umol)")  # Midazolam row of Table 6 column 'K_i' = +Inf (not a reversible inhibitor), entered as 0
    kinactMdz <- fixed(0); label("Midazolam (victim) maximal CYP3A4 inactivation rate k_inact (1/h)")  # Midazolam row of Table 6 column 'k_inact'
    rkiiMdz <- fixed(0); label("Midazolam (victim) reciprocal inactivation constant 1/K_I (L/umol)")  # Midazolam row of Table 6 column 'K_I' = +Inf (not a mechanism-based inhibitor), entered as 0
    fimaxMdz <- fixed(1); label("Midazolam (victim) maximal CYP3A4 induction fold FI_max")  # Midazolam row of Table 6 column 'FI_max'
    rec50Mdz <- fixed(0); label("Midazolam (victim) reciprocal induction potency 1/EC_50* (L/umol)")  # Midazolam row of Table 6 column 'EC_50*' (the permeability- and f_u,inc-corrected value) = +Inf (not an inducer), entered as 0
    faMdz <- fixed(1); label("Midazolam (victim) fraction absorbed")  # Midazolam row of Table 7 column 'F_a'
    kaMdz <- fixed(1.16); label("Midazolam (victim) first-order absorption rate constant (1/h)")  # Midazolam row of Table 7 column 'k_a'
    qgMdz <- fixed(15.44); label("Midazolam (victim) hybrid enterocyte-to-portal-vein flow Q_g (L/h)")  # Midazolam row of Table 7 column 'Q_g'
    clintRMdz <- fixed(0.0421787); label("Midazolam (victim) intrinsic renal clearance (L/h)")  # Eq 22 from Table 3 CL_R = 0.09

    # ---- Perpetrator, default ketoconazole ----
    fubPerp <- fixed(0.0136); label("Perpetrator, default ketoconazole blood fraction unbound")  # ketoconazole row of Table 4 column 'f_u^b'
    fuhPerp <- fixed(0.0075); label("Perpetrator, default ketoconazole hepatocyte fraction unbound")  # ketoconazole row of Table 4 column 'f_u^h'
    fugPerp <- fixed(0.0048); label("Perpetrator, default ketoconazole gut-wall fraction unbound")  # ketoconazole row of Table 4 column 'f_u^gw'
    rbpPerp <- fixed(0.7); label("Perpetrator, default ketoconazole blood-to-plasma ratio")  # ketoconazole row of Table 4 column 'R_BP'
    kpRemPerp <- fixed(1.75); label("Perpetrator, default ketoconazole rest-of-body-to-plasma partition coefficient")  # ketoconazole row of Table 5 column 'K_p,RB'
    kpKidPerp <- fixed(1.01); label("Perpetrator, default ketoconazole kidney-to-plasma partition coefficient")  # ketoconazole row of Table 5 column 'K_p,Kidney'
    kpLunPerp <- fixed(0.39); label("Perpetrator, default ketoconazole lung-to-plasma partition coefficient")  # ketoconazole row of Table 5 column 'K_p,Lungs'
    permPerp <- fixed(36000); label("Perpetrator, default ketoconazole hepatocyte membrane permeability (um/h)")  # ketoconazole row of Table 7 column 'P'
    km2Perp <- fixed(1.52); label("Perpetrator, default ketoconazole Michaelis constant of the non-CYP3A4 pathway K_m,2 (umol/L)")  # ketoconazole row of Table 3 column 'K_m'
    rkm1Perp <- fixed(0); label("Perpetrator, default ketoconazole reciprocal CYP3A4 Michaelis constant 1/K_m,1 (L/umol)")  # ketoconazole row of Table 3 'K_m' with f_m,3A4 = 0; not a CYP3A4 substrate, so K_m,1 = +Inf per the Models section note, entered as 0
    cli3a4Perp <- fixed(0); label("Perpetrator, default ketoconazole CYP3A4 intrinsic clearance per unit hepatocyte volume k_cat*E_0/K_m,1 (1/h)")  # Eq 21 line 1 as f_m,3A4 * CL*_m,int / V_h = 0 * 51.46 / 1.392, Table 3 and Results V_h
    vmax2Perp <- fixed(56.192); label("Perpetrator, default ketoconazole hepatic non-CYP3A4 V_max,2 (umol/L/h)")  # Eq 21 line 2 as (1 - 0) * 51.46 * 1.52 / 1.392
    vmax2gPerp <- fixed(0.0205414); label("Perpetrator, default ketoconazole gut-wall non-CYP3A4 V_max,2^g (umol/L/h)")  # Eq 21 line 3 from S2 Table CYP fractions and amounts
    rkiPerp <- fixed(166.667); label("Perpetrator, default ketoconazole reciprocal competitive inhibition constant 1/K_i (L/umol)")  # ketoconazole row of Table 6 column 'K_i'
    kinactPerp <- fixed(0); label("Perpetrator, default ketoconazole maximal CYP3A4 inactivation rate k_inact (1/h)")  # ketoconazole row of Table 6 column 'k_inact'
    rkiiPerp <- fixed(0); label("Perpetrator, default ketoconazole reciprocal inactivation constant 1/K_I (L/umol)")  # ketoconazole row of Table 6 column 'K_I' = +Inf (not a mechanism-based inhibitor), entered as 0
    fimaxPerp <- fixed(1); label("Perpetrator, default ketoconazole maximal CYP3A4 induction fold FI_max")  # ketoconazole row of Table 6 column 'FI_max'
    rec50Perp <- fixed(0); label("Perpetrator, default ketoconazole reciprocal induction potency 1/EC_50* (L/umol)")  # ketoconazole row of Table 6 column 'EC_50*' (the permeability- and f_u,inc-corrected value) = +Inf (not an inducer), entered as 0
    faPerp <- fixed(1); label("Perpetrator, default ketoconazole fraction absorbed")  # ketoconazole row of Table 7 column 'F_a'
    kaPerp <- fixed(1); label("Perpetrator, default ketoconazole first-order absorption rate constant (1/h)")  # ketoconazole row of Table 7 column 'k_a'
    qgPerp <- fixed(23.34); label("Perpetrator, default ketoconazole hybrid enterocyte-to-portal-vein flow Q_g (L/h)")  # ketoconazole row of Table 7 column 'Q_g'
    clintRPerp <- fixed(0); label("Perpetrator, default ketoconazole intrinsic renal clearance (L/h)")  # Eq 22 from Table 3 CL_R = 0

  })

  model({
    # =====================================================================
    # 1. LOBULE MESH (Fig 1C algorithm; S1 Code PathLength.m)
    #    5 sinusoid levels of length 344.8/185.3/92.5/46.0/22.5 um (Table 8), total 691.1 um,
    #    discretised into 20 equal intervals (21 nodes) by the method of lines.
    #    vb<k> / vh<k> are the blood / hepatocyte volumes (L) of node k and
    #    sx<k> its blood-hepatocyte exchange area (dm^2) x 1e-5, so that
    #    sx<k> * P (um/h) is a permeability clearance in L/h. The node values
    #    are the Fig 1C algorithm's profile rescaled to the paper's published
    #    totals: V_b = 283 mL, V_h = 1392 mL, S_ex = 10046 dm^2 (Results).
    #    Because flow is conserved across levels, vb<k> * v(x_k) / dx equals
    #    the whole-liver flow qLiv at every node, so the advection term below
    #    needs no per-node velocity constant.
    vb1 <- 0.031079099; vh1 <- 0.12297867; sx1 <- 0.0099121469  # level 1, x = 0.0 um
    vb2 <- 0.028725302; vh2 <- 0.11790933; sx2 <- 0.0093205385  # level 1, x = 34.6 um
    vb3 <- 0.026371505; vh3 <- 0.11283998; sx3 <- 0.0087289302  # level 1, x = 69.1 um
    vb4 <- 0.024017708; vh4 <- 0.10777063; sx4 <- 0.0081373218  # level 1, x = 103.7 um
    vb5 <- 0.021663911; vh5 <- 0.10270128; sx5 <- 0.0075457134  # level 1, x = 138.2 um
    vb6 <- 0.019310114; vh6 <- 0.097631927; sx6 <- 0.0069541051  # level 1, x = 172.8 um
    vb7 <- 0.016956317; vh7 <- 0.092562577; sx7 <- 0.0063624967  # level 1, x = 207.3 um
    vb8 <- 0.01460252; vh8 <- 0.087493228; sx8 <- 0.0057708883  # level 1, x = 241.9 um
    vb9 <- 0.012248723; vh9 <- 0.082423878; sx9 <- 0.00517928  # level 1, x = 276.5 um
    vb10 <- 0.0098949265; vh10 <- 0.077354529; sx10 <- 0.0045876716  # level 1, x = 311.0 um
    vb11 <- 0.016380817; vh11 <- 0.063346257; sx11 <- 0.0051692095  # level 2, x = 345.6 um
    vb12 <- 0.014024495; vh12 <- 0.058271469; sx12 <- 0.0045769664  # level 2, x = 380.1 um
    vb13 <- 0.011668173; vh13 <- 0.053196681; sx13 <- 0.0039847234  # level 2, x = 414.7 um
    vb14 <- 0.0093118504; vh14 <- 0.048121893; sx14 <- 0.0033924803  # level 2, x = 449.2 um
    vb15 <- 0.0069555282; vh15 <- 0.043047105; sx15 <- 0.0028002372  # level 2, x = 483.8 um
    vb16 <- 0.004599206; vh16 <- 0.037972317; sx16 <- 0.0022079942  # level 2, x = 518.4 um
    vb17 <- 0.0066745056; vh17 <- 0.028499132; sx17 <- 0.0022069974  # level 3, x = 552.9 um
    vb18 <- 0.0043080171; vh18 <- 0.023402449; sx18 <- 0.0016121992  # level 3, x = 587.5 um
    vb19 <- 0.0019415286; vh19 <- 0.018305766; sx19 <- 0.0010174009  # level 3, x = 622.0 um
    vb20 <- 0.0017910977; vh20 <- 0.011106204; sx20 <- 0.00072187813  # level 4, x = 656.6 um
    vb21 <- 0.00047465227; vh21 <- 0.0050647101; sx21 <- 0.000270821  # level 5, x = 691.1 um

    # =====================================================================
    # 2. SYSTEMIC BLOOD AND TISSUE CONCENTRATIONS (umol/L of blood/tissue)
    cArtM <- a_arterial / vArt
    cVenM <- a_venous / vVen
    cKidM <- a_kidney / vKid
    cLunM <- a_lung / vLun
    cRemM <- a_remainder / vRem
    cPorM <- a_portal / vPor
    cGutM <- a_gut / vGut
    cArtP <- a_arterial_perpetrator / vArt
    cVenP <- a_venous_perpetrator / vVen
    cKidP <- a_kidney_perpetrator / vKid
    cLunP <- a_lung_perpetrator / vLun
    cRemP <- a_remainder_perpetrator / vRem
    cPorP <- a_portal_perpetrator / vPor
    cGutP <- a_gut_perpetrator / vGut

    # Lobule inlet: flow-weighted mix of hepatic artery and portal vein (Eq 6)
    cInM <- (qHa * cArtM + qPor * cPorM) / (qHa + qPor)
    cInP <- (qHa * cArtP + qPor * cPorP) / (qHa + qPor)

    # =====================================================================
    # 3. LOBULE NODE ALGEBRA
    #    u* = unbound hepatocyte concentration; bet = sum over drugs of
    #    u * (1/K_m,1 + 1/K_i + 1/K_I) (Eq 12 denominator); ebar = free
    #    normalised enzyme E_bar (Eq 12); ind / mbi = the induction and
    #    mechanism-based-inactivation terms of Eq 12; j* = the permeation
    #    flux of Eq 4 / Eq 10; r* = the metabolic rate of Eq 10.
    cbM1 <- sinusoid_slab1 / vb1
    chM1 <- hepatocyte_slab1 / vh1
    cbP1 <- sinusoid_perpetrator_slab1 / vb1
    chP1 <- hepatocyte_perpetrator_slab1 / vh1
    uM1 <- fuhMdz * chM1
    uP1 <- fuhPerp * chP1
    bet1 <- uM1 * (rkm1Mdz + rkiMdz + rkiiMdz) + uP1 * (rkm1Perp + rkiPerp + rkiiPerp)
    ebar1 <- enzyme_liver_slab1 / (1 + bet1)
    ind1 <- (fimaxMdz - 1) * uM1 * rec50Mdz / (1 + uM1 * rec50Mdz) + (fimaxPerp - 1) * uP1 * rec50Perp / (1 + uP1 * rec50Perp)
    mbi1 <- (kinactMdz * rkiiMdz * uM1 + kinactPerp * rkiiPerp * uP1) / (kdegLiv * (1 + bet1))
    jM1 <- sx1 * permMdz * (fubMdz * cbM1 - uM1)
    jP1 <- sx1 * permPerp * (fubPerp * cbP1 - uP1)
    rM1 <- vh1 * ((cli3a4Mdz + kinactMdz * rkiiMdz * e0Liv) * ebar1 * uM1 + vmax2Mdz * uM1 / (km2Mdz + uM1))
    rP1 <- vh1 * ((cli3a4Perp + kinactPerp * rkiiPerp * e0Liv) * ebar1 * uP1 + vmax2Perp * uP1 / (km2Perp + uP1))
    cbM2 <- sinusoid_slab2 / vb2
    chM2 <- hepatocyte_slab2 / vh2
    cbP2 <- sinusoid_perpetrator_slab2 / vb2
    chP2 <- hepatocyte_perpetrator_slab2 / vh2
    uM2 <- fuhMdz * chM2
    uP2 <- fuhPerp * chP2
    bet2 <- uM2 * (rkm1Mdz + rkiMdz + rkiiMdz) + uP2 * (rkm1Perp + rkiPerp + rkiiPerp)
    ebar2 <- enzyme_liver_slab2 / (1 + bet2)
    ind2 <- (fimaxMdz - 1) * uM2 * rec50Mdz / (1 + uM2 * rec50Mdz) + (fimaxPerp - 1) * uP2 * rec50Perp / (1 + uP2 * rec50Perp)
    mbi2 <- (kinactMdz * rkiiMdz * uM2 + kinactPerp * rkiiPerp * uP2) / (kdegLiv * (1 + bet2))
    jM2 <- sx2 * permMdz * (fubMdz * cbM2 - uM2)
    jP2 <- sx2 * permPerp * (fubPerp * cbP2 - uP2)
    rM2 <- vh2 * ((cli3a4Mdz + kinactMdz * rkiiMdz * e0Liv) * ebar2 * uM2 + vmax2Mdz * uM2 / (km2Mdz + uM2))
    rP2 <- vh2 * ((cli3a4Perp + kinactPerp * rkiiPerp * e0Liv) * ebar2 * uP2 + vmax2Perp * uP2 / (km2Perp + uP2))
    cbM3 <- sinusoid_slab3 / vb3
    chM3 <- hepatocyte_slab3 / vh3
    cbP3 <- sinusoid_perpetrator_slab3 / vb3
    chP3 <- hepatocyte_perpetrator_slab3 / vh3
    uM3 <- fuhMdz * chM3
    uP3 <- fuhPerp * chP3
    bet3 <- uM3 * (rkm1Mdz + rkiMdz + rkiiMdz) + uP3 * (rkm1Perp + rkiPerp + rkiiPerp)
    ebar3 <- enzyme_liver_slab3 / (1 + bet3)
    ind3 <- (fimaxMdz - 1) * uM3 * rec50Mdz / (1 + uM3 * rec50Mdz) + (fimaxPerp - 1) * uP3 * rec50Perp / (1 + uP3 * rec50Perp)
    mbi3 <- (kinactMdz * rkiiMdz * uM3 + kinactPerp * rkiiPerp * uP3) / (kdegLiv * (1 + bet3))
    jM3 <- sx3 * permMdz * (fubMdz * cbM3 - uM3)
    jP3 <- sx3 * permPerp * (fubPerp * cbP3 - uP3)
    rM3 <- vh3 * ((cli3a4Mdz + kinactMdz * rkiiMdz * e0Liv) * ebar3 * uM3 + vmax2Mdz * uM3 / (km2Mdz + uM3))
    rP3 <- vh3 * ((cli3a4Perp + kinactPerp * rkiiPerp * e0Liv) * ebar3 * uP3 + vmax2Perp * uP3 / (km2Perp + uP3))
    cbM4 <- sinusoid_slab4 / vb4
    chM4 <- hepatocyte_slab4 / vh4
    cbP4 <- sinusoid_perpetrator_slab4 / vb4
    chP4 <- hepatocyte_perpetrator_slab4 / vh4
    uM4 <- fuhMdz * chM4
    uP4 <- fuhPerp * chP4
    bet4 <- uM4 * (rkm1Mdz + rkiMdz + rkiiMdz) + uP4 * (rkm1Perp + rkiPerp + rkiiPerp)
    ebar4 <- enzyme_liver_slab4 / (1 + bet4)
    ind4 <- (fimaxMdz - 1) * uM4 * rec50Mdz / (1 + uM4 * rec50Mdz) + (fimaxPerp - 1) * uP4 * rec50Perp / (1 + uP4 * rec50Perp)
    mbi4 <- (kinactMdz * rkiiMdz * uM4 + kinactPerp * rkiiPerp * uP4) / (kdegLiv * (1 + bet4))
    jM4 <- sx4 * permMdz * (fubMdz * cbM4 - uM4)
    jP4 <- sx4 * permPerp * (fubPerp * cbP4 - uP4)
    rM4 <- vh4 * ((cli3a4Mdz + kinactMdz * rkiiMdz * e0Liv) * ebar4 * uM4 + vmax2Mdz * uM4 / (km2Mdz + uM4))
    rP4 <- vh4 * ((cli3a4Perp + kinactPerp * rkiiPerp * e0Liv) * ebar4 * uP4 + vmax2Perp * uP4 / (km2Perp + uP4))
    cbM5 <- sinusoid_slab5 / vb5
    chM5 <- hepatocyte_slab5 / vh5
    cbP5 <- sinusoid_perpetrator_slab5 / vb5
    chP5 <- hepatocyte_perpetrator_slab5 / vh5
    uM5 <- fuhMdz * chM5
    uP5 <- fuhPerp * chP5
    bet5 <- uM5 * (rkm1Mdz + rkiMdz + rkiiMdz) + uP5 * (rkm1Perp + rkiPerp + rkiiPerp)
    ebar5 <- enzyme_liver_slab5 / (1 + bet5)
    ind5 <- (fimaxMdz - 1) * uM5 * rec50Mdz / (1 + uM5 * rec50Mdz) + (fimaxPerp - 1) * uP5 * rec50Perp / (1 + uP5 * rec50Perp)
    mbi5 <- (kinactMdz * rkiiMdz * uM5 + kinactPerp * rkiiPerp * uP5) / (kdegLiv * (1 + bet5))
    jM5 <- sx5 * permMdz * (fubMdz * cbM5 - uM5)
    jP5 <- sx5 * permPerp * (fubPerp * cbP5 - uP5)
    rM5 <- vh5 * ((cli3a4Mdz + kinactMdz * rkiiMdz * e0Liv) * ebar5 * uM5 + vmax2Mdz * uM5 / (km2Mdz + uM5))
    rP5 <- vh5 * ((cli3a4Perp + kinactPerp * rkiiPerp * e0Liv) * ebar5 * uP5 + vmax2Perp * uP5 / (km2Perp + uP5))
    cbM6 <- sinusoid_slab6 / vb6
    chM6 <- hepatocyte_slab6 / vh6
    cbP6 <- sinusoid_perpetrator_slab6 / vb6
    chP6 <- hepatocyte_perpetrator_slab6 / vh6
    uM6 <- fuhMdz * chM6
    uP6 <- fuhPerp * chP6
    bet6 <- uM6 * (rkm1Mdz + rkiMdz + rkiiMdz) + uP6 * (rkm1Perp + rkiPerp + rkiiPerp)
    ebar6 <- enzyme_liver_slab6 / (1 + bet6)
    ind6 <- (fimaxMdz - 1) * uM6 * rec50Mdz / (1 + uM6 * rec50Mdz) + (fimaxPerp - 1) * uP6 * rec50Perp / (1 + uP6 * rec50Perp)
    mbi6 <- (kinactMdz * rkiiMdz * uM6 + kinactPerp * rkiiPerp * uP6) / (kdegLiv * (1 + bet6))
    jM6 <- sx6 * permMdz * (fubMdz * cbM6 - uM6)
    jP6 <- sx6 * permPerp * (fubPerp * cbP6 - uP6)
    rM6 <- vh6 * ((cli3a4Mdz + kinactMdz * rkiiMdz * e0Liv) * ebar6 * uM6 + vmax2Mdz * uM6 / (km2Mdz + uM6))
    rP6 <- vh6 * ((cli3a4Perp + kinactPerp * rkiiPerp * e0Liv) * ebar6 * uP6 + vmax2Perp * uP6 / (km2Perp + uP6))
    cbM7 <- sinusoid_slab7 / vb7
    chM7 <- hepatocyte_slab7 / vh7
    cbP7 <- sinusoid_perpetrator_slab7 / vb7
    chP7 <- hepatocyte_perpetrator_slab7 / vh7
    uM7 <- fuhMdz * chM7
    uP7 <- fuhPerp * chP7
    bet7 <- uM7 * (rkm1Mdz + rkiMdz + rkiiMdz) + uP7 * (rkm1Perp + rkiPerp + rkiiPerp)
    ebar7 <- enzyme_liver_slab7 / (1 + bet7)
    ind7 <- (fimaxMdz - 1) * uM7 * rec50Mdz / (1 + uM7 * rec50Mdz) + (fimaxPerp - 1) * uP7 * rec50Perp / (1 + uP7 * rec50Perp)
    mbi7 <- (kinactMdz * rkiiMdz * uM7 + kinactPerp * rkiiPerp * uP7) / (kdegLiv * (1 + bet7))
    jM7 <- sx7 * permMdz * (fubMdz * cbM7 - uM7)
    jP7 <- sx7 * permPerp * (fubPerp * cbP7 - uP7)
    rM7 <- vh7 * ((cli3a4Mdz + kinactMdz * rkiiMdz * e0Liv) * ebar7 * uM7 + vmax2Mdz * uM7 / (km2Mdz + uM7))
    rP7 <- vh7 * ((cli3a4Perp + kinactPerp * rkiiPerp * e0Liv) * ebar7 * uP7 + vmax2Perp * uP7 / (km2Perp + uP7))
    cbM8 <- sinusoid_slab8 / vb8
    chM8 <- hepatocyte_slab8 / vh8
    cbP8 <- sinusoid_perpetrator_slab8 / vb8
    chP8 <- hepatocyte_perpetrator_slab8 / vh8
    uM8 <- fuhMdz * chM8
    uP8 <- fuhPerp * chP8
    bet8 <- uM8 * (rkm1Mdz + rkiMdz + rkiiMdz) + uP8 * (rkm1Perp + rkiPerp + rkiiPerp)
    ebar8 <- enzyme_liver_slab8 / (1 + bet8)
    ind8 <- (fimaxMdz - 1) * uM8 * rec50Mdz / (1 + uM8 * rec50Mdz) + (fimaxPerp - 1) * uP8 * rec50Perp / (1 + uP8 * rec50Perp)
    mbi8 <- (kinactMdz * rkiiMdz * uM8 + kinactPerp * rkiiPerp * uP8) / (kdegLiv * (1 + bet8))
    jM8 <- sx8 * permMdz * (fubMdz * cbM8 - uM8)
    jP8 <- sx8 * permPerp * (fubPerp * cbP8 - uP8)
    rM8 <- vh8 * ((cli3a4Mdz + kinactMdz * rkiiMdz * e0Liv) * ebar8 * uM8 + vmax2Mdz * uM8 / (km2Mdz + uM8))
    rP8 <- vh8 * ((cli3a4Perp + kinactPerp * rkiiPerp * e0Liv) * ebar8 * uP8 + vmax2Perp * uP8 / (km2Perp + uP8))
    cbM9 <- sinusoid_slab9 / vb9
    chM9 <- hepatocyte_slab9 / vh9
    cbP9 <- sinusoid_perpetrator_slab9 / vb9
    chP9 <- hepatocyte_perpetrator_slab9 / vh9
    uM9 <- fuhMdz * chM9
    uP9 <- fuhPerp * chP9
    bet9 <- uM9 * (rkm1Mdz + rkiMdz + rkiiMdz) + uP9 * (rkm1Perp + rkiPerp + rkiiPerp)
    ebar9 <- enzyme_liver_slab9 / (1 + bet9)
    ind9 <- (fimaxMdz - 1) * uM9 * rec50Mdz / (1 + uM9 * rec50Mdz) + (fimaxPerp - 1) * uP9 * rec50Perp / (1 + uP9 * rec50Perp)
    mbi9 <- (kinactMdz * rkiiMdz * uM9 + kinactPerp * rkiiPerp * uP9) / (kdegLiv * (1 + bet9))
    jM9 <- sx9 * permMdz * (fubMdz * cbM9 - uM9)
    jP9 <- sx9 * permPerp * (fubPerp * cbP9 - uP9)
    rM9 <- vh9 * ((cli3a4Mdz + kinactMdz * rkiiMdz * e0Liv) * ebar9 * uM9 + vmax2Mdz * uM9 / (km2Mdz + uM9))
    rP9 <- vh9 * ((cli3a4Perp + kinactPerp * rkiiPerp * e0Liv) * ebar9 * uP9 + vmax2Perp * uP9 / (km2Perp + uP9))
    cbM10 <- sinusoid_slab10 / vb10
    chM10 <- hepatocyte_slab10 / vh10
    cbP10 <- sinusoid_perpetrator_slab10 / vb10
    chP10 <- hepatocyte_perpetrator_slab10 / vh10
    uM10 <- fuhMdz * chM10
    uP10 <- fuhPerp * chP10
    bet10 <- uM10 * (rkm1Mdz + rkiMdz + rkiiMdz) + uP10 * (rkm1Perp + rkiPerp + rkiiPerp)
    ebar10 <- enzyme_liver_slab10 / (1 + bet10)
    ind10 <- (fimaxMdz - 1) * uM10 * rec50Mdz / (1 + uM10 * rec50Mdz) + (fimaxPerp - 1) * uP10 * rec50Perp / (1 + uP10 * rec50Perp)
    mbi10 <- (kinactMdz * rkiiMdz * uM10 + kinactPerp * rkiiPerp * uP10) / (kdegLiv * (1 + bet10))
    jM10 <- sx10 * permMdz * (fubMdz * cbM10 - uM10)
    jP10 <- sx10 * permPerp * (fubPerp * cbP10 - uP10)
    rM10 <- vh10 * ((cli3a4Mdz + kinactMdz * rkiiMdz * e0Liv) * ebar10 * uM10 + vmax2Mdz * uM10 / (km2Mdz + uM10))
    rP10 <- vh10 * ((cli3a4Perp + kinactPerp * rkiiPerp * e0Liv) * ebar10 * uP10 + vmax2Perp * uP10 / (km2Perp + uP10))
    cbM11 <- sinusoid_slab11 / vb11
    chM11 <- hepatocyte_slab11 / vh11
    cbP11 <- sinusoid_perpetrator_slab11 / vb11
    chP11 <- hepatocyte_perpetrator_slab11 / vh11
    uM11 <- fuhMdz * chM11
    uP11 <- fuhPerp * chP11
    bet11 <- uM11 * (rkm1Mdz + rkiMdz + rkiiMdz) + uP11 * (rkm1Perp + rkiPerp + rkiiPerp)
    ebar11 <- enzyme_liver_slab11 / (1 + bet11)
    ind11 <- (fimaxMdz - 1) * uM11 * rec50Mdz / (1 + uM11 * rec50Mdz) + (fimaxPerp - 1) * uP11 * rec50Perp / (1 + uP11 * rec50Perp)
    mbi11 <- (kinactMdz * rkiiMdz * uM11 + kinactPerp * rkiiPerp * uP11) / (kdegLiv * (1 + bet11))
    jM11 <- sx11 * permMdz * (fubMdz * cbM11 - uM11)
    jP11 <- sx11 * permPerp * (fubPerp * cbP11 - uP11)
    rM11 <- vh11 * ((cli3a4Mdz + kinactMdz * rkiiMdz * e0Liv) * ebar11 * uM11 + vmax2Mdz * uM11 / (km2Mdz + uM11))
    rP11 <- vh11 * ((cli3a4Perp + kinactPerp * rkiiPerp * e0Liv) * ebar11 * uP11 + vmax2Perp * uP11 / (km2Perp + uP11))
    cbM12 <- sinusoid_slab12 / vb12
    chM12 <- hepatocyte_slab12 / vh12
    cbP12 <- sinusoid_perpetrator_slab12 / vb12
    chP12 <- hepatocyte_perpetrator_slab12 / vh12
    uM12 <- fuhMdz * chM12
    uP12 <- fuhPerp * chP12
    bet12 <- uM12 * (rkm1Mdz + rkiMdz + rkiiMdz) + uP12 * (rkm1Perp + rkiPerp + rkiiPerp)
    ebar12 <- enzyme_liver_slab12 / (1 + bet12)
    ind12 <- (fimaxMdz - 1) * uM12 * rec50Mdz / (1 + uM12 * rec50Mdz) + (fimaxPerp - 1) * uP12 * rec50Perp / (1 + uP12 * rec50Perp)
    mbi12 <- (kinactMdz * rkiiMdz * uM12 + kinactPerp * rkiiPerp * uP12) / (kdegLiv * (1 + bet12))
    jM12 <- sx12 * permMdz * (fubMdz * cbM12 - uM12)
    jP12 <- sx12 * permPerp * (fubPerp * cbP12 - uP12)
    rM12 <- vh12 * ((cli3a4Mdz + kinactMdz * rkiiMdz * e0Liv) * ebar12 * uM12 + vmax2Mdz * uM12 / (km2Mdz + uM12))
    rP12 <- vh12 * ((cli3a4Perp + kinactPerp * rkiiPerp * e0Liv) * ebar12 * uP12 + vmax2Perp * uP12 / (km2Perp + uP12))
    cbM13 <- sinusoid_slab13 / vb13
    chM13 <- hepatocyte_slab13 / vh13
    cbP13 <- sinusoid_perpetrator_slab13 / vb13
    chP13 <- hepatocyte_perpetrator_slab13 / vh13
    uM13 <- fuhMdz * chM13
    uP13 <- fuhPerp * chP13
    bet13 <- uM13 * (rkm1Mdz + rkiMdz + rkiiMdz) + uP13 * (rkm1Perp + rkiPerp + rkiiPerp)
    ebar13 <- enzyme_liver_slab13 / (1 + bet13)
    ind13 <- (fimaxMdz - 1) * uM13 * rec50Mdz / (1 + uM13 * rec50Mdz) + (fimaxPerp - 1) * uP13 * rec50Perp / (1 + uP13 * rec50Perp)
    mbi13 <- (kinactMdz * rkiiMdz * uM13 + kinactPerp * rkiiPerp * uP13) / (kdegLiv * (1 + bet13))
    jM13 <- sx13 * permMdz * (fubMdz * cbM13 - uM13)
    jP13 <- sx13 * permPerp * (fubPerp * cbP13 - uP13)
    rM13 <- vh13 * ((cli3a4Mdz + kinactMdz * rkiiMdz * e0Liv) * ebar13 * uM13 + vmax2Mdz * uM13 / (km2Mdz + uM13))
    rP13 <- vh13 * ((cli3a4Perp + kinactPerp * rkiiPerp * e0Liv) * ebar13 * uP13 + vmax2Perp * uP13 / (km2Perp + uP13))
    cbM14 <- sinusoid_slab14 / vb14
    chM14 <- hepatocyte_slab14 / vh14
    cbP14 <- sinusoid_perpetrator_slab14 / vb14
    chP14 <- hepatocyte_perpetrator_slab14 / vh14
    uM14 <- fuhMdz * chM14
    uP14 <- fuhPerp * chP14
    bet14 <- uM14 * (rkm1Mdz + rkiMdz + rkiiMdz) + uP14 * (rkm1Perp + rkiPerp + rkiiPerp)
    ebar14 <- enzyme_liver_slab14 / (1 + bet14)
    ind14 <- (fimaxMdz - 1) * uM14 * rec50Mdz / (1 + uM14 * rec50Mdz) + (fimaxPerp - 1) * uP14 * rec50Perp / (1 + uP14 * rec50Perp)
    mbi14 <- (kinactMdz * rkiiMdz * uM14 + kinactPerp * rkiiPerp * uP14) / (kdegLiv * (1 + bet14))
    jM14 <- sx14 * permMdz * (fubMdz * cbM14 - uM14)
    jP14 <- sx14 * permPerp * (fubPerp * cbP14 - uP14)
    rM14 <- vh14 * ((cli3a4Mdz + kinactMdz * rkiiMdz * e0Liv) * ebar14 * uM14 + vmax2Mdz * uM14 / (km2Mdz + uM14))
    rP14 <- vh14 * ((cli3a4Perp + kinactPerp * rkiiPerp * e0Liv) * ebar14 * uP14 + vmax2Perp * uP14 / (km2Perp + uP14))
    cbM15 <- sinusoid_slab15 / vb15
    chM15 <- hepatocyte_slab15 / vh15
    cbP15 <- sinusoid_perpetrator_slab15 / vb15
    chP15 <- hepatocyte_perpetrator_slab15 / vh15
    uM15 <- fuhMdz * chM15
    uP15 <- fuhPerp * chP15
    bet15 <- uM15 * (rkm1Mdz + rkiMdz + rkiiMdz) + uP15 * (rkm1Perp + rkiPerp + rkiiPerp)
    ebar15 <- enzyme_liver_slab15 / (1 + bet15)
    ind15 <- (fimaxMdz - 1) * uM15 * rec50Mdz / (1 + uM15 * rec50Mdz) + (fimaxPerp - 1) * uP15 * rec50Perp / (1 + uP15 * rec50Perp)
    mbi15 <- (kinactMdz * rkiiMdz * uM15 + kinactPerp * rkiiPerp * uP15) / (kdegLiv * (1 + bet15))
    jM15 <- sx15 * permMdz * (fubMdz * cbM15 - uM15)
    jP15 <- sx15 * permPerp * (fubPerp * cbP15 - uP15)
    rM15 <- vh15 * ((cli3a4Mdz + kinactMdz * rkiiMdz * e0Liv) * ebar15 * uM15 + vmax2Mdz * uM15 / (km2Mdz + uM15))
    rP15 <- vh15 * ((cli3a4Perp + kinactPerp * rkiiPerp * e0Liv) * ebar15 * uP15 + vmax2Perp * uP15 / (km2Perp + uP15))
    cbM16 <- sinusoid_slab16 / vb16
    chM16 <- hepatocyte_slab16 / vh16
    cbP16 <- sinusoid_perpetrator_slab16 / vb16
    chP16 <- hepatocyte_perpetrator_slab16 / vh16
    uM16 <- fuhMdz * chM16
    uP16 <- fuhPerp * chP16
    bet16 <- uM16 * (rkm1Mdz + rkiMdz + rkiiMdz) + uP16 * (rkm1Perp + rkiPerp + rkiiPerp)
    ebar16 <- enzyme_liver_slab16 / (1 + bet16)
    ind16 <- (fimaxMdz - 1) * uM16 * rec50Mdz / (1 + uM16 * rec50Mdz) + (fimaxPerp - 1) * uP16 * rec50Perp / (1 + uP16 * rec50Perp)
    mbi16 <- (kinactMdz * rkiiMdz * uM16 + kinactPerp * rkiiPerp * uP16) / (kdegLiv * (1 + bet16))
    jM16 <- sx16 * permMdz * (fubMdz * cbM16 - uM16)
    jP16 <- sx16 * permPerp * (fubPerp * cbP16 - uP16)
    rM16 <- vh16 * ((cli3a4Mdz + kinactMdz * rkiiMdz * e0Liv) * ebar16 * uM16 + vmax2Mdz * uM16 / (km2Mdz + uM16))
    rP16 <- vh16 * ((cli3a4Perp + kinactPerp * rkiiPerp * e0Liv) * ebar16 * uP16 + vmax2Perp * uP16 / (km2Perp + uP16))
    cbM17 <- sinusoid_slab17 / vb17
    chM17 <- hepatocyte_slab17 / vh17
    cbP17 <- sinusoid_perpetrator_slab17 / vb17
    chP17 <- hepatocyte_perpetrator_slab17 / vh17
    uM17 <- fuhMdz * chM17
    uP17 <- fuhPerp * chP17
    bet17 <- uM17 * (rkm1Mdz + rkiMdz + rkiiMdz) + uP17 * (rkm1Perp + rkiPerp + rkiiPerp)
    ebar17 <- enzyme_liver_slab17 / (1 + bet17)
    ind17 <- (fimaxMdz - 1) * uM17 * rec50Mdz / (1 + uM17 * rec50Mdz) + (fimaxPerp - 1) * uP17 * rec50Perp / (1 + uP17 * rec50Perp)
    mbi17 <- (kinactMdz * rkiiMdz * uM17 + kinactPerp * rkiiPerp * uP17) / (kdegLiv * (1 + bet17))
    jM17 <- sx17 * permMdz * (fubMdz * cbM17 - uM17)
    jP17 <- sx17 * permPerp * (fubPerp * cbP17 - uP17)
    rM17 <- vh17 * ((cli3a4Mdz + kinactMdz * rkiiMdz * e0Liv) * ebar17 * uM17 + vmax2Mdz * uM17 / (km2Mdz + uM17))
    rP17 <- vh17 * ((cli3a4Perp + kinactPerp * rkiiPerp * e0Liv) * ebar17 * uP17 + vmax2Perp * uP17 / (km2Perp + uP17))
    cbM18 <- sinusoid_slab18 / vb18
    chM18 <- hepatocyte_slab18 / vh18
    cbP18 <- sinusoid_perpetrator_slab18 / vb18
    chP18 <- hepatocyte_perpetrator_slab18 / vh18
    uM18 <- fuhMdz * chM18
    uP18 <- fuhPerp * chP18
    bet18 <- uM18 * (rkm1Mdz + rkiMdz + rkiiMdz) + uP18 * (rkm1Perp + rkiPerp + rkiiPerp)
    ebar18 <- enzyme_liver_slab18 / (1 + bet18)
    ind18 <- (fimaxMdz - 1) * uM18 * rec50Mdz / (1 + uM18 * rec50Mdz) + (fimaxPerp - 1) * uP18 * rec50Perp / (1 + uP18 * rec50Perp)
    mbi18 <- (kinactMdz * rkiiMdz * uM18 + kinactPerp * rkiiPerp * uP18) / (kdegLiv * (1 + bet18))
    jM18 <- sx18 * permMdz * (fubMdz * cbM18 - uM18)
    jP18 <- sx18 * permPerp * (fubPerp * cbP18 - uP18)
    rM18 <- vh18 * ((cli3a4Mdz + kinactMdz * rkiiMdz * e0Liv) * ebar18 * uM18 + vmax2Mdz * uM18 / (km2Mdz + uM18))
    rP18 <- vh18 * ((cli3a4Perp + kinactPerp * rkiiPerp * e0Liv) * ebar18 * uP18 + vmax2Perp * uP18 / (km2Perp + uP18))
    cbM19 <- sinusoid_slab19 / vb19
    chM19 <- hepatocyte_slab19 / vh19
    cbP19 <- sinusoid_perpetrator_slab19 / vb19
    chP19 <- hepatocyte_perpetrator_slab19 / vh19
    uM19 <- fuhMdz * chM19
    uP19 <- fuhPerp * chP19
    bet19 <- uM19 * (rkm1Mdz + rkiMdz + rkiiMdz) + uP19 * (rkm1Perp + rkiPerp + rkiiPerp)
    ebar19 <- enzyme_liver_slab19 / (1 + bet19)
    ind19 <- (fimaxMdz - 1) * uM19 * rec50Mdz / (1 + uM19 * rec50Mdz) + (fimaxPerp - 1) * uP19 * rec50Perp / (1 + uP19 * rec50Perp)
    mbi19 <- (kinactMdz * rkiiMdz * uM19 + kinactPerp * rkiiPerp * uP19) / (kdegLiv * (1 + bet19))
    jM19 <- sx19 * permMdz * (fubMdz * cbM19 - uM19)
    jP19 <- sx19 * permPerp * (fubPerp * cbP19 - uP19)
    rM19 <- vh19 * ((cli3a4Mdz + kinactMdz * rkiiMdz * e0Liv) * ebar19 * uM19 + vmax2Mdz * uM19 / (km2Mdz + uM19))
    rP19 <- vh19 * ((cli3a4Perp + kinactPerp * rkiiPerp * e0Liv) * ebar19 * uP19 + vmax2Perp * uP19 / (km2Perp + uP19))
    cbM20 <- sinusoid_slab20 / vb20
    chM20 <- hepatocyte_slab20 / vh20
    cbP20 <- sinusoid_perpetrator_slab20 / vb20
    chP20 <- hepatocyte_perpetrator_slab20 / vh20
    uM20 <- fuhMdz * chM20
    uP20 <- fuhPerp * chP20
    bet20 <- uM20 * (rkm1Mdz + rkiMdz + rkiiMdz) + uP20 * (rkm1Perp + rkiPerp + rkiiPerp)
    ebar20 <- enzyme_liver_slab20 / (1 + bet20)
    ind20 <- (fimaxMdz - 1) * uM20 * rec50Mdz / (1 + uM20 * rec50Mdz) + (fimaxPerp - 1) * uP20 * rec50Perp / (1 + uP20 * rec50Perp)
    mbi20 <- (kinactMdz * rkiiMdz * uM20 + kinactPerp * rkiiPerp * uP20) / (kdegLiv * (1 + bet20))
    jM20 <- sx20 * permMdz * (fubMdz * cbM20 - uM20)
    jP20 <- sx20 * permPerp * (fubPerp * cbP20 - uP20)
    rM20 <- vh20 * ((cli3a4Mdz + kinactMdz * rkiiMdz * e0Liv) * ebar20 * uM20 + vmax2Mdz * uM20 / (km2Mdz + uM20))
    rP20 <- vh20 * ((cli3a4Perp + kinactPerp * rkiiPerp * e0Liv) * ebar20 * uP20 + vmax2Perp * uP20 / (km2Perp + uP20))
    cbM21 <- sinusoid_slab21 / vb21
    chM21 <- hepatocyte_slab21 / vh21
    cbP21 <- sinusoid_perpetrator_slab21 / vb21
    chP21 <- hepatocyte_perpetrator_slab21 / vh21
    uM21 <- fuhMdz * chM21
    uP21 <- fuhPerp * chP21
    bet21 <- uM21 * (rkm1Mdz + rkiMdz + rkiiMdz) + uP21 * (rkm1Perp + rkiPerp + rkiiPerp)
    ebar21 <- enzyme_liver_slab21 / (1 + bet21)
    ind21 <- (fimaxMdz - 1) * uM21 * rec50Mdz / (1 + uM21 * rec50Mdz) + (fimaxPerp - 1) * uP21 * rec50Perp / (1 + uP21 * rec50Perp)
    mbi21 <- (kinactMdz * rkiiMdz * uM21 + kinactPerp * rkiiPerp * uP21) / (kdegLiv * (1 + bet21))
    jM21 <- sx21 * permMdz * (fubMdz * cbM21 - uM21)
    jP21 <- sx21 * permPerp * (fubPerp * cbP21 - uP21)
    rM21 <- vh21 * ((cli3a4Mdz + kinactMdz * rkiiMdz * e0Liv) * ebar21 * uM21 + vmax2Mdz * uM21 / (km2Mdz + uM21))
    rP21 <- vh21 * ((cli3a4Perp + kinactPerp * rkiiPerp * e0Liv) * ebar21 * uP21 + vmax2Perp * uP21 / (km2Perp + uP21))

    # =====================================================================
    # 4. GUT WALL (Eq 18) -- one well-stirred enterocyte pool per drug
    ugM <- fugMdz * cGutM
    ugP <- fugPerp * cGutP
    betg <- ugM * (rkm1Mdz + rkiMdz + rkiiMdz) + ugP * (rkm1Perp + rkiPerp + rkiiPerp)
    ebarg <- enzyme_gut / (1 + betg)
    indg <- (fimaxMdz - 1) * ugM * rec50Mdz / (1 + ugM * rec50Mdz) + (fimaxPerp - 1) * ugP * rec50Perp / (1 + ugP * rec50Perp)
    mbig <- (kinactMdz * rkiiMdz * ugM + kinactPerp * rkiiPerp * ugP) / (kdegGut * (1 + betg))
    rgM <- vGut * ((cli3a4Mdz * e0Gut / e0Liv + kinactMdz * rkiiMdz * e0Gut) * ebarg * ugM + vmax2gMdz * ugM / (km2Mdz + ugM))
    rgP <- vGut * ((cli3a4Perp * e0Gut / e0Liv + kinactPerp * rkiiPerp * e0Gut) * ebarg * ugP + vmax2gPerp * ugP / (km2Perp + ugP))

    # =====================================================================
    # 5. ODE SYSTEM
    # 5a. Absorption -- the analytical dose sum of Eq 18 written as a depot
    #     state, which is algebraically identical and lets rxode2 handle the
    #     dosing records.
    d/dt(depot) <- -kaMdz * depot
    d/dt(depot_perpetrator) <- -kaPerp * depot_perpetrator
    f(depot) <- faMdz
    f(depot_perpetrator) <- faPerp

    # 5b. Gut wall and portal vein (Eq 18, Eq 19)
    d/dt(a_gut) <- kaMdz * depot - rgM - qgMdz * ugM
    d/dt(a_gut_perpetrator) <- kaPerp * depot_perpetrator - rgP - qgPerp * ugP
    d/dt(a_portal) <- qPor * (cArtM - cPorM) + qgMdz * ugM
    d/dt(a_portal_perpetrator) <- qPor * (cArtP - cPorP) + qgPerp * ugP
    d/dt(enzyme_gut) <- kdegGut * (1 + indg - enzyme_gut * (1 + mbig))

    # 5c. Liver lobule (Eq 4 blood, Eq 10 hepatocytes, Eq 12 enzyme).
    #     Advection is the upwind scheme of S1 Code f1_Mat; node 1 takes the
    #     inlet boundary condition C_0(t) of Eq 6.
    d/dt(sinusoid_slab1) <- qLiv * (cInM - cbM1) - jM1
    d/dt(sinusoid_perpetrator_slab1) <- qLiv * (cInP - cbP1) - jP1
    d/dt(hepatocyte_slab1) <- jM1 - rM1
    d/dt(hepatocyte_perpetrator_slab1) <- jP1 - rP1
    d/dt(enzyme_liver_slab1) <- kdegLiv * (1 + ind1 - enzyme_liver_slab1 * (1 + mbi1))
    d/dt(sinusoid_slab2) <- qLiv * (cbM1 - cbM2) - jM2
    d/dt(sinusoid_perpetrator_slab2) <- qLiv * (cbP1 - cbP2) - jP2
    d/dt(hepatocyte_slab2) <- jM2 - rM2
    d/dt(hepatocyte_perpetrator_slab2) <- jP2 - rP2
    d/dt(enzyme_liver_slab2) <- kdegLiv * (1 + ind2 - enzyme_liver_slab2 * (1 + mbi2))
    d/dt(sinusoid_slab3) <- qLiv * (cbM2 - cbM3) - jM3
    d/dt(sinusoid_perpetrator_slab3) <- qLiv * (cbP2 - cbP3) - jP3
    d/dt(hepatocyte_slab3) <- jM3 - rM3
    d/dt(hepatocyte_perpetrator_slab3) <- jP3 - rP3
    d/dt(enzyme_liver_slab3) <- kdegLiv * (1 + ind3 - enzyme_liver_slab3 * (1 + mbi3))
    d/dt(sinusoid_slab4) <- qLiv * (cbM3 - cbM4) - jM4
    d/dt(sinusoid_perpetrator_slab4) <- qLiv * (cbP3 - cbP4) - jP4
    d/dt(hepatocyte_slab4) <- jM4 - rM4
    d/dt(hepatocyte_perpetrator_slab4) <- jP4 - rP4
    d/dt(enzyme_liver_slab4) <- kdegLiv * (1 + ind4 - enzyme_liver_slab4 * (1 + mbi4))
    d/dt(sinusoid_slab5) <- qLiv * (cbM4 - cbM5) - jM5
    d/dt(sinusoid_perpetrator_slab5) <- qLiv * (cbP4 - cbP5) - jP5
    d/dt(hepatocyte_slab5) <- jM5 - rM5
    d/dt(hepatocyte_perpetrator_slab5) <- jP5 - rP5
    d/dt(enzyme_liver_slab5) <- kdegLiv * (1 + ind5 - enzyme_liver_slab5 * (1 + mbi5))
    d/dt(sinusoid_slab6) <- qLiv * (cbM5 - cbM6) - jM6
    d/dt(sinusoid_perpetrator_slab6) <- qLiv * (cbP5 - cbP6) - jP6
    d/dt(hepatocyte_slab6) <- jM6 - rM6
    d/dt(hepatocyte_perpetrator_slab6) <- jP6 - rP6
    d/dt(enzyme_liver_slab6) <- kdegLiv * (1 + ind6 - enzyme_liver_slab6 * (1 + mbi6))
    d/dt(sinusoid_slab7) <- qLiv * (cbM6 - cbM7) - jM7
    d/dt(sinusoid_perpetrator_slab7) <- qLiv * (cbP6 - cbP7) - jP7
    d/dt(hepatocyte_slab7) <- jM7 - rM7
    d/dt(hepatocyte_perpetrator_slab7) <- jP7 - rP7
    d/dt(enzyme_liver_slab7) <- kdegLiv * (1 + ind7 - enzyme_liver_slab7 * (1 + mbi7))
    d/dt(sinusoid_slab8) <- qLiv * (cbM7 - cbM8) - jM8
    d/dt(sinusoid_perpetrator_slab8) <- qLiv * (cbP7 - cbP8) - jP8
    d/dt(hepatocyte_slab8) <- jM8 - rM8
    d/dt(hepatocyte_perpetrator_slab8) <- jP8 - rP8
    d/dt(enzyme_liver_slab8) <- kdegLiv * (1 + ind8 - enzyme_liver_slab8 * (1 + mbi8))
    d/dt(sinusoid_slab9) <- qLiv * (cbM8 - cbM9) - jM9
    d/dt(sinusoid_perpetrator_slab9) <- qLiv * (cbP8 - cbP9) - jP9
    d/dt(hepatocyte_slab9) <- jM9 - rM9
    d/dt(hepatocyte_perpetrator_slab9) <- jP9 - rP9
    d/dt(enzyme_liver_slab9) <- kdegLiv * (1 + ind9 - enzyme_liver_slab9 * (1 + mbi9))
    d/dt(sinusoid_slab10) <- qLiv * (cbM9 - cbM10) - jM10
    d/dt(sinusoid_perpetrator_slab10) <- qLiv * (cbP9 - cbP10) - jP10
    d/dt(hepatocyte_slab10) <- jM10 - rM10
    d/dt(hepatocyte_perpetrator_slab10) <- jP10 - rP10
    d/dt(enzyme_liver_slab10) <- kdegLiv * (1 + ind10 - enzyme_liver_slab10 * (1 + mbi10))
    d/dt(sinusoid_slab11) <- qLiv * (cbM10 - cbM11) - jM11
    d/dt(sinusoid_perpetrator_slab11) <- qLiv * (cbP10 - cbP11) - jP11
    d/dt(hepatocyte_slab11) <- jM11 - rM11
    d/dt(hepatocyte_perpetrator_slab11) <- jP11 - rP11
    d/dt(enzyme_liver_slab11) <- kdegLiv * (1 + ind11 - enzyme_liver_slab11 * (1 + mbi11))
    d/dt(sinusoid_slab12) <- qLiv * (cbM11 - cbM12) - jM12
    d/dt(sinusoid_perpetrator_slab12) <- qLiv * (cbP11 - cbP12) - jP12
    d/dt(hepatocyte_slab12) <- jM12 - rM12
    d/dt(hepatocyte_perpetrator_slab12) <- jP12 - rP12
    d/dt(enzyme_liver_slab12) <- kdegLiv * (1 + ind12 - enzyme_liver_slab12 * (1 + mbi12))
    d/dt(sinusoid_slab13) <- qLiv * (cbM12 - cbM13) - jM13
    d/dt(sinusoid_perpetrator_slab13) <- qLiv * (cbP12 - cbP13) - jP13
    d/dt(hepatocyte_slab13) <- jM13 - rM13
    d/dt(hepatocyte_perpetrator_slab13) <- jP13 - rP13
    d/dt(enzyme_liver_slab13) <- kdegLiv * (1 + ind13 - enzyme_liver_slab13 * (1 + mbi13))
    d/dt(sinusoid_slab14) <- qLiv * (cbM13 - cbM14) - jM14
    d/dt(sinusoid_perpetrator_slab14) <- qLiv * (cbP13 - cbP14) - jP14
    d/dt(hepatocyte_slab14) <- jM14 - rM14
    d/dt(hepatocyte_perpetrator_slab14) <- jP14 - rP14
    d/dt(enzyme_liver_slab14) <- kdegLiv * (1 + ind14 - enzyme_liver_slab14 * (1 + mbi14))
    d/dt(sinusoid_slab15) <- qLiv * (cbM14 - cbM15) - jM15
    d/dt(sinusoid_perpetrator_slab15) <- qLiv * (cbP14 - cbP15) - jP15
    d/dt(hepatocyte_slab15) <- jM15 - rM15
    d/dt(hepatocyte_perpetrator_slab15) <- jP15 - rP15
    d/dt(enzyme_liver_slab15) <- kdegLiv * (1 + ind15 - enzyme_liver_slab15 * (1 + mbi15))
    d/dt(sinusoid_slab16) <- qLiv * (cbM15 - cbM16) - jM16
    d/dt(sinusoid_perpetrator_slab16) <- qLiv * (cbP15 - cbP16) - jP16
    d/dt(hepatocyte_slab16) <- jM16 - rM16
    d/dt(hepatocyte_perpetrator_slab16) <- jP16 - rP16
    d/dt(enzyme_liver_slab16) <- kdegLiv * (1 + ind16 - enzyme_liver_slab16 * (1 + mbi16))
    d/dt(sinusoid_slab17) <- qLiv * (cbM16 - cbM17) - jM17
    d/dt(sinusoid_perpetrator_slab17) <- qLiv * (cbP16 - cbP17) - jP17
    d/dt(hepatocyte_slab17) <- jM17 - rM17
    d/dt(hepatocyte_perpetrator_slab17) <- jP17 - rP17
    d/dt(enzyme_liver_slab17) <- kdegLiv * (1 + ind17 - enzyme_liver_slab17 * (1 + mbi17))
    d/dt(sinusoid_slab18) <- qLiv * (cbM17 - cbM18) - jM18
    d/dt(sinusoid_perpetrator_slab18) <- qLiv * (cbP17 - cbP18) - jP18
    d/dt(hepatocyte_slab18) <- jM18 - rM18
    d/dt(hepatocyte_perpetrator_slab18) <- jP18 - rP18
    d/dt(enzyme_liver_slab18) <- kdegLiv * (1 + ind18 - enzyme_liver_slab18 * (1 + mbi18))
    d/dt(sinusoid_slab19) <- qLiv * (cbM18 - cbM19) - jM19
    d/dt(sinusoid_perpetrator_slab19) <- qLiv * (cbP18 - cbP19) - jP19
    d/dt(hepatocyte_slab19) <- jM19 - rM19
    d/dt(hepatocyte_perpetrator_slab19) <- jP19 - rP19
    d/dt(enzyme_liver_slab19) <- kdegLiv * (1 + ind19 - enzyme_liver_slab19 * (1 + mbi19))
    d/dt(sinusoid_slab20) <- qLiv * (cbM19 - cbM20) - jM20
    d/dt(sinusoid_perpetrator_slab20) <- qLiv * (cbP19 - cbP20) - jP20
    d/dt(hepatocyte_slab20) <- jM20 - rM20
    d/dt(hepatocyte_perpetrator_slab20) <- jP20 - rP20
    d/dt(enzyme_liver_slab20) <- kdegLiv * (1 + ind20 - enzyme_liver_slab20 * (1 + mbi20))
    d/dt(sinusoid_slab21) <- qLiv * (cbM20 - cbM21) - jM21
    d/dt(sinusoid_perpetrator_slab21) <- qLiv * (cbP20 - cbP21) - jP21
    d/dt(hepatocyte_slab21) <- jM21 - rM21
    d/dt(hepatocyte_perpetrator_slab21) <- jP21 - rP21
    d/dt(enzyme_liver_slab21) <- kdegLiv * (1 + ind21 - enzyme_liver_slab21 * (1 + mbi21))

    # 5d. Systemic compartments (Eq 13-17). The liver outlet is the blood
    #     concentration at the last lobule node (node 21).
    d/dt(a_arterial) <- qTot * (cLunM * rbpMdz / kpLunMdz - cArtM)
    d/dt(a_arterial_perpetrator) <- qTot * (cLunP * rbpPerp / kpLunPerp - cArtP)
    d/dt(a_venous) <- qLiv * cbM21 + qKid * cKidM * rbpMdz / kpKidMdz + qRem * cRemM * rbpMdz / kpRemMdz - qTot * cVenM
    d/dt(a_venous_perpetrator) <- qLiv * cbP21 + qKid * cKidP * rbpPerp / kpKidPerp + qRem * cRemP * rbpPerp / kpRemPerp - qTot * cVenP
    d/dt(a_kidney) <- qKid * (cArtM - cKidM * rbpMdz / kpKidMdz) - clintRMdz * cKidM
    d/dt(a_kidney_perpetrator) <- qKid * (cArtP - cKidP * rbpPerp / kpKidPerp) - clintRPerp * cKidP
    d/dt(a_lung) <- qTot * (cVenM - cLunM * rbpMdz / kpLunMdz)
    d/dt(a_lung_perpetrator) <- qTot * (cVenP - cLunP * rbpPerp / kpLunPerp)
    d/dt(a_remainder) <- qRem * (cArtM - cRemM * rbpMdz / kpRemMdz)
    d/dt(a_remainder_perpetrator) <- qRem * (cArtP - cRemP * rbpPerp / kpRemPerp)

    # 5e. Enzyme states start at their basal level (Eq 12 normalisation)
    enzyme_gut(0) <- 1
    enzyme_liver_slab1(0) <- 1
    enzyme_liver_slab2(0) <- 1
    enzyme_liver_slab3(0) <- 1
    enzyme_liver_slab4(0) <- 1
    enzyme_liver_slab5(0) <- 1
    enzyme_liver_slab6(0) <- 1
    enzyme_liver_slab7(0) <- 1
    enzyme_liver_slab8(0) <- 1
    enzyme_liver_slab9(0) <- 1
    enzyme_liver_slab10(0) <- 1
    enzyme_liver_slab11(0) <- 1
    enzyme_liver_slab12(0) <- 1
    enzyme_liver_slab13(0) <- 1
    enzyme_liver_slab14(0) <- 1
    enzyme_liver_slab15(0) <- 1
    enzyme_liver_slab16(0) <- 1
    enzyme_liver_slab17(0) <- 1
    enzyme_liver_slab18(0) <- 1
    enzyme_liver_slab19(0) <- 1
    enzyme_liver_slab20(0) <- 1
    enzyme_liver_slab21(0) <- 1

    # =====================================================================
    # 6. OBSERVATION -- midazolam PLASMA concentration (umol/L). The
    #    systemic states hold BLOOD concentrations, so divide by R_BP.
    Cpp <- cVenP / rbpPerp
    Cc <- cVenM / rbpMdz
  })
}
