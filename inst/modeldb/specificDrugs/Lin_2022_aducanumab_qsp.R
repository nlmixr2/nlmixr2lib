Lin_2022_aducanumab_qsp <- function() {
  description <- paste(
    "QSP. Amyloid-beta (A-beta) pathway and aducanumab mechanism-of-action",
    "model in early Alzheimer's disease. Thirty-six mass-action ODE states in",
    "three physiologic compartments (plasma 3 L, CSF 0.139 L, brain",
    "interstitial fluid 0.261 L) plus a peripheral aducanumab compartment:",
    "APP synthesis and sequential beta-secretase (BACE) and gamma-secretase",
    "cleavage to A-beta monomer in plasma and brain ISF, shedding of soluble",
    "BACE, monomer-oligomer exchange in every compartment, oligomer-plaque",
    "exchange in brain ISF, inter-compartment transport of monomer, oligomer,",
    "soluble BACE, drug and soluble drug-A-beta complexes, two-compartment",
    "aducanumab disposition, drug binding to monomer, oligomer and plaque,",
    "and FcR-mediated antibody-dependent cellular phagocytosis (ADCP) that",
    "clears drug-oligomer and drug-plaque complexes in brain ISF. States are",
    "amounts in nmol; second-order rates divide by the volume of the",
    "compartment in which the reaction occurs. The primary PD output is the",
    "percent change in total brain plaque, which the paper equates with the",
    "percent change in amyloid-PET composite SUVR above a cutoff of 1.0.",
    "Deterministic: no between-subject variability or residual error is",
    "reported. Nominal (final) parameter set; the paper's alternative",
    "faster-plaque-turnover set (Table S3) is reproduced in the vignette.",
    "The system is stiff and badly scaled (states span 1e-3 to 1.5e5 nmol):",
    "solve with tight tolerances, e.g. rxSolve(..., atol = 1e-14, rtol =",
    "1e-10, maxsteps = 5e6); rxode2's defaults fail.",
    sep = " "
  )
  reference <- paste(
    "Lin L, Hua F, Salinas C, Young C, Bussiere T, Apgar JF, Burke JM,",
    "Kandadi Muralidharan K, Rajagovindan R, Nestorov I. Quantitative",
    "systems pharmacology model for Alzheimer's disease to predict the",
    "effect of aducanumab on brain amyloid. CPT Pharmacometrics Syst",
    "Pharmacol. 2022;11(3):362-372. doi:10.1002/psp4.12759.",
    "Reaction network, compartment volumes and parameter values from",
    "Supplementary Model Code (S2, KroneckerBio model file); parameter",
    "provenance and units from Table S2; pretreatment steady state from the",
    "Supplementary_Initial_Condition sheet (S3); equations from S4.",
    sep = " "
  )
  vignette <- "Lin_2022_aducanumab_qsp"

  # The published reaction network (Supplementary Model Code) resolves each
  # A-beta pathway species and each drug complex in each physiologic
  # compartment. None of these maps onto a canonical nlmixr2lib compartment
  # role except free drug in plasma (`central`) and the aducanumab
  # peripheral compartment (`peripheral1`). The `_plasma` / `_csf` / `_bisf`
  # suffix follows the paper's compartment names; `bisf` is brain ISF.
  paper_specific_compartments <- c(
    "app_plasma",
    "bace_plasma",
    "app_bace_plasma",
    "ctfb_plasma",
    "gamma_plasma",
    "ctfb_gamma_plasma",
    "baces_plasma",
    "abeta_plasma",
    "aolig_plasma",
    "mab_abeta_plasma",
    "mab_aolig_plasma",
    "mab_csf",
    "abeta_csf",
    "aolig_csf",
    "mab_abeta_csf",
    "mab_aolig_csf",
    "baces_csf",
    "mab_bisf",
    "app_bisf",
    "bace_bisf",
    "app_bace_bisf",
    "ctfb_bisf",
    "gamma_bisf",
    "ctfb_gamma_bisf",
    "baces_bisf",
    "abeta_bisf",
    "aolig_bisf",
    "aplaq_bisf",
    "mab_abeta_bisf",
    "mab_aolig_bisf",
    "mab_aplaq_bisf",
    "fcr_bisf",
    "mab_aolig_fcr_bisf",
    "mab_aplaq_fcr_bisf"
  )

  units <- list(
    time = "day",
    dosing = "mg",
    concentration = paste(
      "Cc and Ccsf are total (free + A-beta-bound) aducanumab in ug/mL;",
      "abeta_plasma_total is total plasma A-beta in pg/mL; plaque_bisf and",
      "aolig_free_bisf are brain-ISF concentrations in nM; pct_plaque and",
      "pct_aolig_free are percent change from the pretreatment steady state.",
      sep = " "
    )
  )

  # Every state holds an amount in nmol (KroneckerBio mass-action
  # convention of the deposited model file: states are amounts, a reaction
  # is attached to a compartment, and bimolecular rates are k * x1 * x2 / V).
  compartmentData <- list(
    central = list(analyte = "aducanumab (free)", units = "nmol", specimen = "plasma", verified = TRUE),
    app_plasma = list(
      analyte = "amyloid precursor protein (APP)",
      units = "nmol",
      specimen = "plasma",
      verified = TRUE
    ),
    bace_plasma = list(
      analyte = "membrane beta-secretase (BACE)",
      units = "nmol",
      specimen = "plasma",
      verified = TRUE
    ),
    app_bace_plasma = list(
      analyte = "APP-BACE enzyme-substrate complex",
      units = "nmol",
      specimen = "plasma",
      verified = TRUE
    ),
    ctfb_plasma = list(
      analyte = "C-terminal fragment beta (CTF-beta)",
      units = "nmol",
      specimen = "plasma",
      verified = TRUE
    ),
    gamma_plasma = list(analyte = "gamma-secretase", units = "nmol", specimen = "plasma", verified = TRUE),
    ctfb_gamma_plasma = list(
      analyte = "CTF-beta-gamma-secretase enzyme-substrate complex",
      units = "nmol",
      specimen = "plasma",
      verified = TRUE
    ),
    baces_plasma = list(
      analyte = "soluble beta-secretase (sBACE)",
      units = "nmol",
      specimen = "plasma",
      verified = TRUE
    ),
    abeta_plasma = list(
      analyte = "amyloid-beta monomer (A-beta 1-40 + 1-42)",
      units = "nmol",
      specimen = "plasma",
      verified = TRUE
    ),
    aolig_plasma = list(
      analyte = "soluble amyloid-beta oligomer",
      units = "nmol",
      specimen = "plasma",
      verified = TRUE
    ),
    mab_abeta_plasma = list(
      analyte = "aducanumab-A-beta monomer complex",
      units = "nmol",
      specimen = "plasma",
      verified = TRUE
    ),
    mab_aolig_plasma = list(
      analyte = "aducanumab-A-beta oligomer complex",
      units = "nmol",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(analyte = "aducanumab (free)", units = "nmol", specimen = "tissue", verified = TRUE),
    mab_csf = list(analyte = "aducanumab (free)", units = "nmol", specimen = "CSF", verified = TRUE),
    abeta_csf = list(
      analyte = "amyloid-beta monomer (A-beta 1-40 + 1-42)",
      units = "nmol",
      specimen = "CSF",
      verified = TRUE
    ),
    aolig_csf = list(analyte = "soluble amyloid-beta oligomer", units = "nmol", specimen = "CSF", verified = TRUE),
    mab_abeta_csf = list(
      analyte = "aducanumab-A-beta monomer complex",
      units = "nmol",
      specimen = "CSF",
      verified = TRUE
    ),
    mab_aolig_csf = list(
      analyte = "aducanumab-A-beta oligomer complex",
      units = "nmol",
      specimen = "CSF",
      verified = TRUE
    ),
    baces_csf = list(analyte = "soluble beta-secretase (sBACE)", units = "nmol", specimen = "CSF", verified = TRUE),
    mab_bisf = list(analyte = "aducanumab (free)", units = "nmol", specimen = "brain ISF", verified = TRUE),
    app_bisf = list(
      analyte = "amyloid precursor protein (APP)",
      units = "nmol",
      specimen = "brain ISF",
      verified = TRUE
    ),
    bace_bisf = list(
      analyte = "membrane beta-secretase (BACE)",
      units = "nmol",
      specimen = "brain ISF",
      verified = TRUE
    ),
    app_bace_bisf = list(
      analyte = "APP-BACE enzyme-substrate complex",
      units = "nmol",
      specimen = "brain ISF",
      verified = TRUE
    ),
    ctfb_bisf = list(
      analyte = "C-terminal fragment beta (CTF-beta)",
      units = "nmol",
      specimen = "brain ISF",
      verified = TRUE
    ),
    gamma_bisf = list(analyte = "gamma-secretase", units = "nmol", specimen = "brain ISF", verified = TRUE),
    ctfb_gamma_bisf = list(
      analyte = "CTF-beta-gamma-secretase enzyme-substrate complex",
      units = "nmol",
      specimen = "brain ISF",
      verified = TRUE
    ),
    baces_bisf = list(
      analyte = "soluble beta-secretase (sBACE)",
      units = "nmol",
      specimen = "brain ISF",
      verified = TRUE
    ),
    abeta_bisf = list(
      analyte = "amyloid-beta monomer (A-beta 1-40 + 1-42)",
      units = "nmol",
      specimen = "brain ISF",
      verified = TRUE
    ),
    aolig_bisf = list(
      analyte = "soluble amyloid-beta oligomer",
      units = "nmol",
      specimen = "brain ISF",
      verified = TRUE
    ),
    aplaq_bisf = list(
      analyte = "insoluble amyloid-beta plaque",
      units = "nmol",
      specimen = "brain ISF",
      verified = TRUE
    ),
    mab_abeta_bisf = list(
      analyte = "aducanumab-A-beta monomer complex",
      units = "nmol",
      specimen = "brain ISF",
      verified = TRUE
    ),
    mab_aolig_bisf = list(
      analyte = "aducanumab-A-beta oligomer complex",
      units = "nmol",
      specimen = "brain ISF",
      verified = TRUE
    ),
    mab_aplaq_bisf = list(
      analyte = "aducanumab-A-beta plaque complex",
      units = "nmol",
      specimen = "brain ISF",
      verified = TRUE
    ),
    fcr_bisf = list(
      analyte = "Fc-gamma receptor on phagocytic (microglial) cells",
      units = "nmol",
      specimen = "brain ISF",
      verified = TRUE
    ),
    mab_aolig_fcr_bisf = list(
      analyte = "aducanumab-A-beta oligomer-FcR ternary complex",
      units = "nmol",
      specimen = "brain ISF",
      verified = TRUE
    ),
    mab_aplaq_fcr_bisf = list(
      analyte = "aducanumab-A-beta plaque-FcR ternary complex",
      units = "nmol",
      specimen = "brain ISF",
      verified = TRUE
    )
  )

  covariateData <- list()

  covariatesDataExcluded <- list()

  population <- list(
    species = "human",
    n_subjects = NA_integer_,
    n_studies = 3L,
    disease_state = paste(
      "Mild-to-moderate Alzheimer's disease (single-ascending-dose study) and",
      "prodromal or mild Alzheimer's disease with amyloid-PET-positive scans",
      "(phase Ib PRIME multiple-ascending-dose study, including its long-term",
      "extension and dose-titration cohort)"
    ),
    dose_range = paste(
      "Calibration: 0.3-60 mg/kg single i.v. dose (SAD) and 1-10 mg/kg i.v.",
      "q4w for 1 year (MAD placebo-controlled period). Validation: 1-10 mg/kg",
      "q4w to week 110 (LTE) and a 1 -> 3 -> 6 -> 10 mg/kg q4w titration."
    ),
    regions = "Not reported (Biogen clinical program)",
    notes = paste(
      "QSP calibrated to group-mean data, not a population fit. Table S1",
      "lists the calibration and validation data: literature baseline",
      "concentrations of A-beta monomer, oligomer and plaque in plasma, CSF",
      "and brain ISF; SILK A-beta kinetics in CSF (Mawuenyega 2010); SAD",
      "serum PK and plasma A-beta (Ferrero 2016); CSF aducanumab (internal);",
      "1-year SUVR from the PRIME placebo-controlled period (Sevigny 2016);",
      "validation against MAD PK, week-110 LTE SUVR and the titration cohort",
      "(Viglietta 2017). Subject counts and demographics are not reported in",
      "this paper. The model represents a typical 70 kg patient (inferred",
      "from Figure 2a; see the vignette).",
      sep = " "
    )
  )

  ini({
    # =====================================================================
    # Units. Every rate constant is published in 1/s, nmol/s or 1/nM/s
    # (Table S2; Supplementary Model Code). Each is multiplied by 86400 here
    # so the model runs on a time base of DAYS; the published per-second
    # value is kept visible inside log() so the source trace reads directly
    # against the supplement.
    #
    # Values are those of the deposited KroneckerBio model file
    # (Supplementary Model Code, % Parameters), which prints more digits
    # than Table S2. Where the two disagree, the value used here is the one
    # that reproduces the paper's own published outputs; each case is
    # annotated on its line and discussed in the vignette:
    #   ksynthAPP_plasma, ksynthBACE_plasma -- recovered from the published
    #     steady state (S3), which the code's 4-decimal display rounds;
    #   kG2M_bisf    -- code 1.40e-8 (S3 steady state), not Table S2 1.48e-8;
    #   koffma0      -- code 10 /s (Figure 2b), not Table S2 1.00 /s;
    #   koffma1/2    -- Table S2 and text 2.00e-2 /s (Figures 2c, 3a),
    #                   not code 0.0180 /s.
    #
    # fixed() marks values Table S2 attributes to the literature or to an
    # assumption; values Table S2 marks 'Fitted' in the calibration are
    # left unfixed. No uncertainty is reported for any value.
    # =====================================================================

    # ---- Physiologic volumes --------------------------------------------
    lvc <- fixed(log(3)); label("Plasma volume (L)") # Supplementary Model Code, % Compartments: vplasma 3
    lvcsf <- fixed(log(0.139)); label("CSF volume (L)") # Supplementary Model Code, % Compartments: vcsf 0.139
    lvbisf <- fixed(log(0.261)); label("Brain interstitial fluid volume (L)") # Supplementary Model Code, % Compartments: vbisf 0.261

    # ---- APP synthesis and degradation ----------------------------------
    lkclearapp <- log(4.8135e-05 * 86400); label("APP degradation rate constant kclearAPP (1/day)") # Code kclearAPP 4.8135e-05 /s; Table S2 4.81e-5 /s, 'Fitted to SILK data'
    # Code displays 0.0014 at 4 decimals; the S3 steady state (APP_plasma
    # 28.93565 nmol) gives kclearAPP*APP + kcatBACE*APP_BACE = 1.44405e-3.
    lksynthapp_plasma <- fixed(log(1.44405e-03 * 86400)); label("APP synthesis rate in plasma ksynthAPP_plasma (nmol/day)") # Code 0.0014 nmol/s (display-rounded); Table S2 1.40e-3; value back-solved from S3 steady state
    lksynthapp_bisf <- fixed(log(2.5127e-04 * 86400)); label("APP synthesis rate in brain ISF ksynthAPP_bisf (nmol/day)") # Code 2.5127e-04 nmol/s; Table S2 2.51e-4 (refs 1,3)

    # ---- BACE synthesis, shedding and degradation -----------------------
    lkclearbace <- fixed(log(1.2034e-05 * 86400)); label("BACE degradation rate constant kclearBACE (1/day)") # Code 1.2034e-05 /s; Table S2 1.20e-5 (ref 11)
    # Code displays 0.0011; the S3 steady state (BACE_plasma 88.52867 nmol)
    # gives (kclearBACE + kcleave) * BACE = 1.08306e-3.
    lksynthbace_plasma <- log(1.08306e-03 * 86400); label("BACE synthesis rate in plasma ksynthBACE_plasma (nmol/day)") # Code 0.0011 nmol/s (display-rounded); Table S2 1.10e-3 'Fitted'; value back-solved from S3 steady state
    lksynthbace_bisf <- log(1.0365e-04 * 86400); label("BACE synthesis rate in brain ISF ksynthBACE_bisf (nmol/day)") # Code 1.0365e-04 nmol/s; Table S2 1.04e-04 'Fitted'
    lkcleave <- log(2.0e-07 * 86400); label("Cleavage of BACE to soluble BACE kcleave (1/day)") # Code 2.0000e-07 /s; Table S2 2.00e-7 'Fitted to BACEs data'
    lkclearbaces <- log(6.4180e-05 * 86400); label("Soluble BACE clearance rate constant kclearBACEs (1/day)") # Code 6.4180e-05 /s; Table S2 6.42e-5 'Fitted to BACEs data'

    # ---- Gamma-secretase synthesis and degradation ----------------------
    lkcleargamma <- log(1.9254e-04 * 86400); label("Gamma-secretase clearance rate constant kclearGamma (1/day)") # Code 1.9254e-04 /s; Table S2 1.93e-4 'Fitted to A-beta data'
    lksynthgamma_plasma <- log(28.8811 * 86400); label("Gamma-secretase synthesis rate in plasma ksynthGamma_plasma (nmol/day)") # Code 28.8811 nmol/s; Table S2 28.9
    lksynthgamma_bisf <- log(2.5127 * 86400); label("Gamma-secretase synthesis rate in brain ISF ksynthGamma_bisf (nmol/day)") # Code ksynthGamma_bisf 2.5127 nmol/s; Table S2 2.51 (row mislabelled 'ksynthGamma_plasma')

    # ---- Enzyme binding and catalysis -----------------------------------
    lkonpp <- fixed(log(1.0e-03 * 86400)); label("Protein-protein association rate constant konPP (1/nM/day)") # Code konPP 1.0000e-03 /nM/s; Table S2 1.00e-3, typical value (ref 15)
    lkoffbace <- fixed(log(119.9928 * 86400)); label("APP-BACE dissociation rate constant koffBACE (1/day)") # Code 119.9928 /s; Table S2 120 (refs 16-18)
    lkcatbace <- fixed(log(0.0072 * 86400)); label("BACE catalytic rate constant kcatBACE (1/day)") # Code 0.0072 /s; Table S2 0.0072 (ref 19)
    lkoffgamma <- fixed(log(215.9983 * 86400)); label("CTF-beta-gamma-secretase dissociation rate constant koffGamma (1/day)") # Code 215.9983 /s; Table S2 216 (ref 20)
    lkcatgamma <- fixed(log(0.0017 * 86400)); label("Gamma-secretase catalytic rate constant kcatGamma (1/day)") # Code 0.0017 /s; Table S2 0.0017 (ref 20)

    # ---- A-beta clearance -----------------------------------------------
    lkclearabeta_plasma <- fixed(log(9.6270e-05 * 86400)); label("A-beta monomer clearance in plasma kclearAbeta_plasma (1/day)") # Code 9.6270e-05 /s; Table S2 9.63e-5 (ref 14)
    lkclearaolig_plasma <- fixed(log(9.6270e-05 * 86400)); label("A-beta oligomer clearance in plasma kclearAolig_plasma (1/day)") # Code 9.6270e-05 /s; Table S2 9.63e-5
    lkclearabeta_bisf <- fixed(log(1.9254e-05 * 86400)); label("A-beta monomer degradation in brain ISF kclearAbeta_bisf (1/day)") # Code 1.9254e-05 /s; Table S2 1.93e-05 (ref 27)
    lkclearaolig_bisf <- log(2.2040e-08 * 86400); label("A-beta oligomer degradation in brain ISF kclearAolig_bisf (1/day)") # Code 2.2040e-08 /s; Table S2 2.20e-08 'Fitted to A-beta baseline data'
    lkclearaplaq <- fixed(log(4.4080e-09 * 86400)); label("Endogenous plaque degradation in brain ISF kclearAplaq (1/day)") # Code 4.4080e-09 /s; Table S2 4.41e-09 (ref 28)

    # ---- Aggregation ----------------------------------------------------
    lkm2g <- fixed(log(1.4e-05 * 86400)); label("Monomer-to-oligomer rate constant kM2G (1/day)") # Code 1.4000e-05 /s; Table S2 1.40e-5 (refs 21,22)
    lkg2m_plasma <- fixed(log(140000 * 86400)); label("Oligomer-to-monomer rate constant in plasma kG2M_plasma (1/day)") # Code 140000 /s; Table S2 1.40E+05, assumed large
    lkg2m_csf <- log(0.0028 * 86400); label("Oligomer-to-monomer rate constant in CSF kG2M_csf (1/day)") # Code 0.0028 /s; Table S2 2.80e-03 'Fit A-beta baseline data in CSF'
    lkg2m_bisf <- log(1.4e-08 * 86400); label("Oligomer-to-monomer rate constant in brain ISF kG2M_bisf (1/day)") # Code 1.4000e-08 /s (matches S3 steady state); Table S2 prints 1.48e-08
    lkg2p <- log(7.0e-08 * 86400); label("Oligomer-to-plaque rate constant kG2P (1/day)") # Code 7.0000e-08 /s; Table S2 7.00e-08 'Fit A-beta baseline data in brain'
    lkp2g_bisf <- fixed(log(7.0e-11 * 86400)); label("Plaque-to-oligomer rate constant kP2G_bisf (1/day)") # Code 7.0000e-11 /s; Table S2 7.00e-11, assumed

    # ---- Transport plasma <-> CSF (index 1 = plasma, 3 = CSF) -----------
    lk13abeta <- fixed(log(1.7222e-09 * 86400)); label("A-beta monomer transport plasma to CSF k13Abeta (1/day)") # Code 1.7222e-09 /s; Table S2 1.72e-09 (refs 7,23-25)
    lk31abeta <- fixed(log(4.1667e-05 * 86400)); label("A-beta monomer transport CSF to plasma k31Abeta (1/day)") # Code 4.1667e-05 /s; Table S2 4.17e-05
    lk13aolig <- fixed(log(1.7222e-09 * 86400)); label("A-beta oligomer transport plasma to CSF k13Aolig (1/day)") # Code 1.7222e-09 /s; Table S2 1.72e-09
    lk31aolig <- fixed(log(4.1667e-05 * 86400)); label("A-beta oligomer transport CSF to plasma k31Aolig (1/day)") # Code 4.1667e-05 /s; Table S2 4.17e-05
    lk13bace <- fixed(log(1.7222e-09 * 86400)); label("Soluble BACE transport plasma to CSF k13BACE (1/day)") # Code 1.7222e-09 /s; Table S2 1.72e-09
    lk31bace <- fixed(log(4.1667e-05 * 86400)); label("Soluble BACE transport CSF to plasma k31BACE (1/day)") # Code 4.1667e-05 /s; Table S2 4.17e-05
    lk13mab <- fixed(log(1.7222e-09 * 86400)); label("Drug and drug-monomer transport plasma to CSF k13mAb (1/day)") # Code 1.7222e-09 /s; Table S2 1.72e-09 (refs 23-25)
    lk31mab <- fixed(log(4.1667e-05 * 86400)); label("Drug and drug-monomer transport CSF to plasma k31mAb (1/day)") # Code 4.1667e-05 /s; Table S2 4.17e-05
    lk13mix <- fixed(log(1.7222e-09 * 86400)); label("Drug-oligomer transport plasma to CSF k13mix (1/day)") # Code k13mix 1.7222e-09 /s (not tabulated separately in Table S2)
    lk31mix <- fixed(log(4.1667e-05 * 86400)); label("Drug-oligomer transport CSF to plasma k31mix (1/day)") # Code k31mix 4.1667e-05 /s (not tabulated separately in Table S2)

    # ---- Transport plasma <-> brain ISF (index 4 = brain ISF) -----------
    lk14abeta <- fixed(log(1.4811e-04 * 86400)); label("A-beta monomer transport plasma to brain ISF k14Abeta (1/day)") # Code 1.4811e-04 /s; Table S2 1.48e-04, scaled from mice (ref 26)
    lk41abeta <- log(1.4811e-05 * 86400); label("A-beta monomer transport brain ISF to plasma k41Abeta (1/day)") # Code 1.4811e-05 /s; Table S2 1.48e-05 'Fitted to A-beta baseline data'
    lk14aolig <- log(1.4811e-06 * 86400); label("A-beta oligomer transport plasma to brain ISF k14Aolig (1/day)") # Code 1.4811e-06 /s; Table S2 1.48e-06
    lk41aolig <- log(1.4811e-08 * 86400); label("A-beta oligomer transport brain ISF to plasma k41Aolig (1/day)") # Code 1.4811e-08 /s; Table S2 1.48e-08
    lk14bace <- log(2.0056e-08 * 86400); label("Soluble BACE transport plasma to brain ISF k14BACE (1/day)") # Code 2.0056e-08 /s; Table S2 2.01e-08 'Fitted to BACEs data in CSF'
    lk41bace <- log(8.0226e-05 * 86400); label("Soluble BACE transport brain ISF to plasma k41BACE (1/day)") # Code 8.0226e-05 /s; Table S2 8.02e-05
    lk14mab <- log(1.6045e-06 * 86400); label("Drug and drug-monomer transport plasma to brain ISF k14mAb (1/day)") # Code 1.6045e-06 /s; Table S2 1.60e-06 'Fitted to CSF aducanumab data'
    lk41mab <- log(0.0032 * 86400); label("Drug and drug-monomer transport brain ISF to plasma k41mAb (1/day)") # Code 0.0032 /s; Table S2 3.20e-03 'Fitted to CSF aducanumab data'
    lk14mix <- fixed(log(1.6045e-06 * 86400)); label("Drug-oligomer transport plasma to brain ISF k14mix (1/day)") # Code 1.6045e-06 /s; Table S2 1.60e-06, assumed
    lk41mix <- fixed(log(1.4811e-08 * 86400)); label("Drug-oligomer transport brain ISF to plasma k41mix (1/day)") # Code 1.4811e-08 /s; Table S2 1.48e-08, assumed equal to oligomer

    # ---- Transport brain ISF -> CSF (one-way) ---------------------------
    lk43abeta <- log(1.5509e-05 * 86400); label("A-beta monomer transport brain ISF to CSF k43Abeta (1/day)") # Code 1.5509e-05 /s; Table S2 1.55e-05 'Fitted to A-beta baseline and SILK data'
    lk43aolig <- log(2.3264e-08 * 86400); label("A-beta oligomer transport brain ISF to CSF k43Aolig (1/day)") # Code 2.3264e-08 /s; Table S2 2.33e-08 'Fitted to A-beta baseline data'
    lk43bace <- fixed(log(1.5509e-05 * 86400)); label("Soluble BACE transport brain ISF to CSF k43BACE (1/day)") # Code 1.5509e-05 /s; Table S2 1.55e-05, assumed equal to monomer
    lk43mab <- log(1.5509e-05 * 86400); label("Drug and drug-monomer transport brain ISF to CSF k43mAb (1/day)") # Code 1.5509e-05 /s; Table S2 1.55e-05 'Fitted to CSF aducanumab data'
    lk43mix <- fixed(log(2.3264e-08 * 86400)); label("Drug-oligomer transport brain ISF to CSF k43mix (1/day)") # Code 2.3264e-08 /s; Table S2 2.33e-08, assumed equal to oligomer

    # ---- Aducanumab disposition -----------------------------------------
    lkclearmab <- log(1.4586e-06 * 86400); label("First-order elimination of drug and drug complexes in plasma kclearmAb (1/day)") # Code 1.4586e-06 /s; Table S2 1.46e-06 'Fitted to aducanumab PK data'
    lk12 <- log(2.5e-06 * 86400); label("Drug distribution plasma to peripheral k12mAb (1/day)") # Code k12mAb 2.5000e-06 /s; Table S2 2.50e-6
    lk21 <- log(1.0e-06 * 86400); label("Drug distribution peripheral to plasma k21mAb (1/day)") # Code k21mAb 1.0000e-06 /s; Table S2 1.00e-06

    # ---- Drug binding to A-beta -----------------------------------------
    lkoffma0 <- log(10 * 86400); label("Drug-A-beta monomer dissociation rate constant koffma0 (1/day)") # Code koffma0 10 /s; Table S2 prints 1.00 /s -- 10 reproduces Figure 2b
    lkoffma1 <- fixed(log(2.0e-02 * 86400)); label("Drug-A-beta oligomer dissociation rate constant koffma1 (1/day)") # Table S2 2.00e-02 /s, assumed equal to plaque; code 0.0180 -- 0.02 reproduces Figures 2c, 3a
    lkonpd <- fixed(log(1.0e-03 * 86400)); label("Drug-plaque association rate constant konPD (1/nM/day)") # Code konPD 1.0000e-03 /nM/s; Table S2 1.00e-3 (ref 15)
    lkoffma2 <- log(2.0e-02 * 86400); label("Drug-plaque dissociation rate constant koffma2 (1/day)") # Table S2 2.00e-2 /s 'Fitted to SUVR data'; text Kd 20 nM; code 0.0180 -- 0.02 reproduces Figures 2c, 3a

    # ---- FcR-mediated phagocytosis (ADCP) -------------------------------
    lkclearfcr <- fixed(log(1.9254e-04 * 86400)); label("FcR degradation rate constant kclearFcR (1/day)") # Code 1.9254e-04 /s; Table S2 1.93e-04, typical receptor turnover (ref 29)
    lksynthfcr_bisf <- fixed(log(5.0253e-05 * 86400)); label("FcR synthesis rate in brain ISF ksynthFcR_bisf (nmol/day)") # Code 5.0253e-05 nmol/s; Table S2 5.03e-05, assumed high
    lkonpf <- fixed(log(1.0e-03 * 86400)); label("Drug-A-beta complex to FcR association rate constant konPF (1/nM/day)") # Code konPF 1.0000e-03 /nM/s; Table S2 1.00e-3 (ref 15)
    lkoffpf <- fixed(log(10 * 86400)); label("Drug-A-beta complex to FcR dissociation rate constant koffPF (1/day)") # Code koffPF 10 /s; Table S2 10.0, assumed (linear ADCP)
    lkcatadcp <- log(0.0036 * 86400); label("ADCP catalytic rate constant kcatADCP (1/day)") # Code 0.0036 /s; Table S2 3.60e-3 'Fitted to SUVR data'

    # ---- Unit conversions -----------------------------------------------
    # Neither molar mass is printed in the paper. Both were back-solved by
    # the maintainers from the paper's own figures; see the vignette.
    mw_adu <- fixed(150000); label("Aducanumab molar mass used for mg <-> nmol conversion (g/mol)") # Back-solved: Figure 3a curves at 70 kg (RMS 0.3 %-points vs 1.0 at 145912 g/mol)
    mw_abeta <- fixed(4330); label("A-beta molar mass used for nM -> pg/mL conversion (g/mol)") # Back-solved: Figure 2b baseline 500 pg/mL = S3 plasma A-beta 0.1155 nM (A-beta 1-40, 4330 g/mol)
  })

  model({
    # =====================================================================
    # 1. Parameters on the natural scale (per day)
    # =====================================================================
    vc <- exp(lvc)
    vcsf <- exp(lvcsf)
    vbisf <- exp(lvbisf)

    kclearapp <- exp(lkclearapp)
    ksynthapp_plasma <- exp(lksynthapp_plasma)
    ksynthapp_bisf <- exp(lksynthapp_bisf)
    kclearbace <- exp(lkclearbace)
    ksynthbace_plasma <- exp(lksynthbace_plasma)
    ksynthbace_bisf <- exp(lksynthbace_bisf)
    kcleave <- exp(lkcleave)
    kclearbaces <- exp(lkclearbaces)
    kcleargamma <- exp(lkcleargamma)
    ksynthgamma_plasma <- exp(lksynthgamma_plasma)
    ksynthgamma_bisf <- exp(lksynthgamma_bisf)
    konpp <- exp(lkonpp)
    koffbace <- exp(lkoffbace)
    kcatbace <- exp(lkcatbace)
    koffgamma <- exp(lkoffgamma)
    kcatgamma <- exp(lkcatgamma)
    kclearabeta_plasma <- exp(lkclearabeta_plasma)
    kclearaolig_plasma <- exp(lkclearaolig_plasma)
    kclearabeta_bisf <- exp(lkclearabeta_bisf)
    kclearaolig_bisf <- exp(lkclearaolig_bisf)
    kclearaplaq <- exp(lkclearaplaq)
    km2g <- exp(lkm2g)
    kg2m_plasma <- exp(lkg2m_plasma)
    kg2m_csf <- exp(lkg2m_csf)
    kg2m_bisf <- exp(lkg2m_bisf)
    kg2p <- exp(lkg2p)
    kp2g_bisf <- exp(lkp2g_bisf)
    k13abeta <- exp(lk13abeta)
    k31abeta <- exp(lk31abeta)
    k13aolig <- exp(lk13aolig)
    k31aolig <- exp(lk31aolig)
    k13bace <- exp(lk13bace)
    k31bace <- exp(lk31bace)
    k13mab <- exp(lk13mab)
    k31mab <- exp(lk31mab)
    k13mix <- exp(lk13mix)
    k31mix <- exp(lk31mix)
    k14abeta <- exp(lk14abeta)
    k41abeta <- exp(lk41abeta)
    k14aolig <- exp(lk14aolig)
    k41aolig <- exp(lk41aolig)
    k14bace <- exp(lk14bace)
    k41bace <- exp(lk41bace)
    k14mab <- exp(lk14mab)
    k41mab <- exp(lk41mab)
    k14mix <- exp(lk14mix)
    k41mix <- exp(lk41mix)
    k43abeta <- exp(lk43abeta)
    k43aolig <- exp(lk43aolig)
    k43bace <- exp(lk43bace)
    k43mab <- exp(lk43mab)
    k43mix <- exp(lk43mix)
    kclearmab <- exp(lkclearmab)
    k12 <- exp(lk12)
    k21 <- exp(lk21)
    koffma0 <- exp(lkoffma0)
    koffma1 <- exp(lkoffma1)
    konpd <- exp(lkonpd)
    koffma2 <- exp(lkoffma2)
    kclearfcr <- exp(lkclearfcr)
    ksynthfcr_bisf <- exp(lksynthfcr_bisf)
    konpf <- exp(lkonpf)
    koffpf <- exp(lkoffpf)
    kcatadcp <- exp(lkcatadcp)

    # =====================================================================
    # 2. Reaction fluxes (nmol/day). Numbering follows the order of
    #    % Reactions in the Supplementary Model Code. Zero-order synthesis
    #    is an amount rate; first-order terms act on amounts; bimolecular
    #    terms are k * x1 * x2 / V with V the compartment of the reaction.
    # =====================================================================
    # -- Plasma -----------------------------------------------------------
    # The code writes APP synthesis as APP_IN -> APP_plasma with an input
    # APP_IN; the input is 1 in every simulation of the paper except the
    # SILK labelling calibration, so it is written here as the constant rate.
    r_app_bace_p <- konpp * app_plasma * bace_plasma / vc - koffbace * app_bace_plasma
    r_ctfb_gamma_p <- konpp * ctfb_plasma * gamma_plasma / vc - koffgamma * ctfb_gamma_plasma
    r_m2g_p <- km2g * abeta_plasma - kg2m_plasma * aolig_plasma
    r_mab_abeta_p <- konpp * central * abeta_plasma / vc - koffma0 * mab_abeta_plasma
    r_mab_aolig_p <- konpp * central * aolig_plasma / vc - koffma1 * mab_aolig_plasma
    # -- Plasma <-> CSF ---------------------------------------------------
    t_mab_pc <- k13mab * central - k31mab * mab_csf
    t_abeta_pc <- k13abeta * abeta_plasma - k31abeta * abeta_csf
    t_aolig_pc <- k13aolig * aolig_plasma - k31aolig * aolig_csf
    t_mab_abeta_pc <- k13mab * mab_abeta_plasma - k31mab * mab_abeta_csf
    t_mab_aolig_pc <- k13mix * mab_aolig_plasma - k31mix * mab_aolig_csf
    t_baces_pc <- k13bace * baces_plasma - k31bace * baces_csf
    # -- Plasma <-> brain ISF ---------------------------------------------
    t_mab_pb <- k14mab * central - k41mab * mab_bisf
    t_abeta_pb <- k14abeta * abeta_plasma - k41abeta * abeta_bisf
    t_aolig_pb <- k14aolig * aolig_plasma - k41aolig * aolig_bisf
    t_mab_abeta_pb <- k14mab * mab_abeta_plasma - k41mab * mab_abeta_bisf
    t_mab_aolig_pb <- k14mix * mab_aolig_plasma - k41mix * mab_aolig_bisf
    t_baces_pb <- k14bace * baces_plasma - k41bace * baces_bisf
    # -- CSF ---------------------------------------------------------------
    r_m2g_c <- km2g * abeta_csf - kg2m_csf * aolig_csf
    r_mab_abeta_c <- konpp * mab_csf * abeta_csf / vcsf - koffma0 * mab_abeta_csf
    r_mab_aolig_c <- konpp * mab_csf * aolig_csf / vcsf - koffma1 * mab_aolig_csf
    # -- Brain ISF ---------------------------------------------------------
    r_app_bace_b <- konpp * app_bisf * bace_bisf / vbisf - koffbace * app_bace_bisf
    r_ctfb_gamma_b <- konpp * ctfb_bisf * gamma_bisf / vbisf - koffgamma * ctfb_gamma_bisf
    r_m2g_b <- km2g * abeta_bisf - kg2m_bisf * aolig_bisf
    r_g2p_b <- kg2p * aolig_bisf - kp2g_bisf * aplaq_bisf
    r_mab_abeta_b <- konpp * mab_bisf * abeta_bisf / vbisf - koffma0 * mab_abeta_bisf
    r_mab_aolig_b <- konpp * mab_bisf * aolig_bisf / vbisf - koffma1 * mab_aolig_bisf
    r_mab_aplaq_b <- konpd * mab_bisf * aplaq_bisf / vbisf - koffma2 * mab_aplaq_bisf
    r_fcr_aolig_b <- konpf * mab_aolig_bisf * fcr_bisf / vbisf - koffpf * mab_aolig_fcr_bisf
    r_fcr_aplaq_b <- konpf * mab_aplaq_bisf * fcr_bisf / vbisf - koffpf * mab_aplaq_fcr_bisf
    adcp_aolig <- kcatadcp * mab_aolig_fcr_bisf
    adcp_aplaq <- kcatadcp * mab_aplaq_fcr_bisf
    # -- Brain ISF -> CSF (one-way) ---------------------------------------
    t_mab_bc <- k43mab * mab_bisf
    t_abeta_bc <- k43abeta * abeta_bisf
    t_aolig_bc <- k43aolig * aolig_bisf
    t_mab_abeta_bc <- k43mab * mab_abeta_bisf
    t_mab_aolig_bc <- k43mix * mab_aolig_bisf
    t_baces_bc <- k43bace * baces_bisf

    # =====================================================================
    # 3. ODE system. Free aducanumab in plasma is `central` (dosing
    #    compartment); the aducanumab peripheral compartment is
    #    `peripheral1`.
    # =====================================================================
    # -- Plasma -----------------------------------------------------------
    d/dt(central) <- -kclearmab * central - r_mab_abeta_p - r_mab_aolig_p -
      k12 * central + k21 * peripheral1 - t_mab_pc - t_mab_pb
    d/dt(app_plasma) <- ksynthapp_plasma - kclearapp * app_plasma - r_app_bace_p
    d/dt(bace_plasma) <- ksynthbace_plasma - kclearbace * bace_plasma - kcleave * bace_plasma -
      r_app_bace_p + kcatbace * app_bace_plasma
    d/dt(app_bace_plasma) <- r_app_bace_p - kcatbace * app_bace_plasma
    d/dt(ctfb_plasma) <- kcatbace * app_bace_plasma - r_ctfb_gamma_p
    d/dt(gamma_plasma) <- ksynthgamma_plasma - kcleargamma * gamma_plasma - r_ctfb_gamma_p +
      kcatgamma * ctfb_gamma_plasma
    d/dt(ctfb_gamma_plasma) <- r_ctfb_gamma_p - kcatgamma * ctfb_gamma_plasma
    d/dt(baces_plasma) <- kcleave * bace_plasma - kclearbaces * baces_plasma - t_baces_pc - t_baces_pb
    d/dt(abeta_plasma) <- kcatgamma * ctfb_gamma_plasma - kclearabeta_plasma * abeta_plasma -
      r_m2g_p - t_abeta_pc - t_abeta_pb - r_mab_abeta_p
    d/dt(aolig_plasma) <- -kclearaolig_plasma * aolig_plasma + r_m2g_p - t_aolig_pc -
      t_aolig_pb - r_mab_aolig_p
    d/dt(mab_abeta_plasma) <- -kclearmab * mab_abeta_plasma + r_mab_abeta_p - t_mab_abeta_pc -
      t_mab_abeta_pb
    d/dt(mab_aolig_plasma) <- -kclearmab * mab_aolig_plasma + r_mab_aolig_p - t_mab_aolig_pc -
      t_mab_aolig_pb
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    # -- CSF ---------------------------------------------------------------
    d/dt(mab_csf) <- t_mab_pc + t_mab_bc - r_mab_abeta_c - r_mab_aolig_c
    d/dt(abeta_csf) <- t_abeta_pc + t_abeta_bc - r_m2g_c - r_mab_abeta_c
    d/dt(aolig_csf) <- t_aolig_pc + t_aolig_bc + r_m2g_c - r_mab_aolig_c
    d/dt(mab_abeta_csf) <- t_mab_abeta_pc + t_mab_abeta_bc + r_mab_abeta_c
    d/dt(mab_aolig_csf) <- t_mab_aolig_pc + t_mab_aolig_bc + r_mab_aolig_c
    d/dt(baces_csf) <- t_baces_pc + t_baces_bc
    # -- Brain ISF ---------------------------------------------------------
    d/dt(mab_bisf) <- t_mab_pb - t_mab_bc - r_mab_abeta_b - r_mab_aolig_b - r_mab_aplaq_b
    d/dt(app_bisf) <- ksynthapp_bisf - kclearapp * app_bisf - r_app_bace_b
    d/dt(bace_bisf) <- ksynthbace_bisf - kclearbace * bace_bisf - kcleave * bace_bisf -
      r_app_bace_b + kcatbace * app_bace_bisf
    d/dt(app_bace_bisf) <- r_app_bace_b - kcatbace * app_bace_bisf
    d/dt(ctfb_bisf) <- kcatbace * app_bace_bisf - r_ctfb_gamma_b
    d/dt(gamma_bisf) <- ksynthgamma_bisf - kcleargamma * gamma_bisf - r_ctfb_gamma_b +
      kcatgamma * ctfb_gamma_bisf
    d/dt(ctfb_gamma_bisf) <- r_ctfb_gamma_b - kcatgamma * ctfb_gamma_bisf
    d/dt(baces_bisf) <- kcleave * bace_bisf - kclearbaces * baces_bisf + t_baces_pb - t_baces_bc
    d/dt(abeta_bisf) <- kcatgamma * ctfb_gamma_bisf - kclearabeta_bisf * abeta_bisf - r_m2g_b +
      t_abeta_pb - t_abeta_bc - r_mab_abeta_b
    d/dt(aolig_bisf) <- -kclearaolig_bisf * aolig_bisf + r_m2g_b - r_g2p_b + t_aolig_pb -
      t_aolig_bc - r_mab_aolig_b
    d/dt(aplaq_bisf) <- r_g2p_b - kclearaplaq * aplaq_bisf - r_mab_aplaq_b
    d/dt(mab_abeta_bisf) <- t_mab_abeta_pb - t_mab_abeta_bc + r_mab_abeta_b
    d/dt(mab_aolig_bisf) <- t_mab_aolig_pb - t_mab_aolig_bc + r_mab_aolig_b - r_fcr_aolig_b
    d/dt(mab_aplaq_bisf) <- r_mab_aplaq_b - r_fcr_aplaq_b
    d/dt(fcr_bisf) <- ksynthfcr_bisf - kclearfcr * fcr_bisf - r_fcr_aolig_b - r_fcr_aplaq_b +
      adcp_aolig + adcp_aplaq
    d/dt(mab_aolig_fcr_bisf) <- r_fcr_aolig_b - adcp_aolig
    d/dt(mab_aplaq_fcr_bisf) <- r_fcr_aplaq_b - adcp_aplaq

    # =====================================================================
    # 4. Initial conditions: pretreatment steady state, from the
    #    Supplementary_Initial_Condition sheet (S3), in nmol. Every drug
    #    and drug-complex state is 0 (S3 lists values below 1e-38, which are
    #    numerical residue). With the parameters above this is a steady
    #    state to within 1e-4 relative over 1000 years of untreated
    #    simulation (vignette). It is the steady state of THIS parameter
    #    set only -- after changing a parameter, run the model untreated to
    #    a new steady state before dosing.
    # =====================================================================
    app_plasma(0) <- 28.93565064
    bace_plasma(0) <- 88.52866711
    app_bace_plasma(0) <- 7.115652e-3
    ctfb_plasma(0) <- 0.130191311
    abeta_plasma(0) <- 0.346458044
    aolig_plasma(0) <- 4.57e-11
    gamma_plasma(0) <- 150000
    ctfb_gamma_plasma(0) <- 3.0136877e-2
    baces_plasma(0) <- 0.291643718
    abeta_csf(0) <- 0.416367423
    aolig_csf(0) <- 2.838285e-3
    baces_csf(0) <- 3.969583e-3
    app_bisf(0) <- 5.017003728
    bace_bisf(0) <- 8.472193442
    app_bace_bisf(0) <- 1.357121e-3
    ctfb_bisf(0) <- 2.4830522e-2
    abeta_bisf(0) <- 0.982018951
    aolig_bisf(0) <- 96.12773625
    aplaq_bisf(0) <- 1502.537143
    gamma_bisf(0) <- 13050
    ctfb_gamma_bisf(0) <- 5.747806e-3
    baces_bisf(0) <- 1.0632428e-2
    fcr_bisf(0) <- 0.261

    # =====================================================================
    # 5. Dosing: doses are given in mg as a 2-hour i.v. infusion into
    #    `central` (Supplementary Model Code, kinfusion = dose in nmol per
    #    2 hours); f() converts mg to nmol.
    # =====================================================================
    f(central) <- 1e6 / mw_adu

    # =====================================================================
    # 6. Outputs (Supplementary Model Code, % Outputs)
    # =====================================================================
    # totalmAbPlasma, total aducanumab in plasma (ug/mL), Figure 2a
    Cc <- (central + mab_abeta_plasma + mab_aolig_plasma) / vc * mw_adu * 1e-6
    # totalmAbCSF, total aducanumab in CSF (ug/mL), Figure S2
    Ccsf <- (mab_csf + mab_abeta_csf + mab_aolig_csf) / vcsf * mw_adu * 1e-6
    # totalAbetaPlasma, total A-beta in plasma (pg/mL), Figure 2b
    abeta_plasma_total <- (abeta_plasma + mab_abeta_plasma + aolig_plasma + mab_aolig_plasma) /
      vc * mw_abeta
    # totalPlaque, total brain plaque (nM), and its percent change from the
    # pretreatment steady state (S3 Aplaq_bisf 1502.537143 nmol, no
    # complexes). The paper equates this with the percent change in
    # composite SUVR above a cutoff of 1.0 (Methods, 'SUVR data
    # processing'); Figures 2c, 3, 4 and 6.
    plaque_bisf <- (aplaq_bisf + mab_aplaq_bisf + mab_aplaq_fcr_bisf) / vbisf
    pct_plaque <- 100 * ((aplaq_bisf + mab_aplaq_bisf + mab_aplaq_fcr_bisf) / 1502.537143 - 1)
    # Free (unbound) soluble oligomer in brain ISF (nM) and its percent
    # change from the S3 steady state (96.12773625 nmol); Figure 5.
    aolig_free_bisf <- aolig_bisf / vbisf
    pct_aolig_free <- 100 * (aolig_bisf / 96.12773625 - 1)

    # Deterministic QSP: the paper reports no between-subject variability
    # and no residual-error model, so none is attached.
  })
}
