Yan_2026_rilzabrutinib_qsp <- function() {
  description <- paste(
    "QSP. BTK target-engagement (BTK occupancy, BO) model for rilzabrutinib",
    "(RIL), a reversible covalent Bruton's tyrosine kinase inhibitor, in",
    "healthy adults. Rilzabrutinib binds reversibly to a BTK pool of",
    "0.1 uM; BO, the bound fraction of BTK, is the paper's efficacy",
    "biomarker (target BOmax > 90% and BOtrough >= 70% at steady state).",
    "The rilzabrutinib concentration is NOT modelled here: Yan 2026",
    "generated it with a whole-body PK-Sim 12.1 PBPK model that is a",
    "platform port (Rodgers and Rowland distribution, PK-Sim database",
    "physiology, no organ volumes, blood flows or partition coefficients",
    "published and no project file deposited), so that layer -- including",
    "its CYP3A4 drug-drug-interaction and hepatic-impairment scenarios -- is",
    "deliberately NOT extracted. Total plasma concentration is supplied per",
    "record as the time-varying covariate CP_RILZABRUTINIB_NGML and the",
    "fraction unbound in plasma (0.025) is applied inside model().",
    "KNOWN GAP: driven by free plasma concentration at the paper's own",
    "predicted exposure (Table 2), the model reproduces the shape of the",
    "published BO time course (rise, 16 h complex half-life, OD versus BID",
    "ordering) but sits about 25 BO points below the paper's predicted and",
    "observed values (Table 3); about 3.8-fold free plasma (effective",
    "unbound fraction about 0.094) is needed to reproduce them. The paper",
    "does not say in which PK-Sim matrix the BTK pool sits, so the",
    "concentration that drove its BO is not published; see the vignette.",
    "Deterministic mechanism model: no IIV and no residual error.",
    sep = " "
  )
  reference <- paste(
    "Yan X, Yan Z, Xiao L, Lu D. (2026). PBPK modeling predicts",
    "rilzabrutinib dose adjustments using BTK occupancy for CYP3A4",
    "interactions and hepatic impairment. Sci Rep 16:20275.",
    "doi:10.1038/s41598-026-50990-5. PMID 42067571, PMCID PMC13323760.",
    "BO equation from Methods Eq. 1 (adapted from Wong 2019, ref 21);",
    "parameter values from Table 1; hepatic-impairment protein-binding",
    "scale factor from Methods Eq. 2 with albumin fractions from",
    "Supplementary Table S3.",
    sep = " "
  )
  vignette <- "Yan_2026_rilzabrutinib"

  units <- list(
    time = "h",
    dosing = "(not applicable; the rilzabrutinib concentration is supplied as the time-varying covariate CP_RILZABRUTINIB_NGML)",
    concentration = "nmol/L"
  )

  # BO was measured in whole blood (Supplementary Table S1), but Yan 2026
  # does not say in which PK-Sim matrix the BTK pool sits, so the specimen
  # is the readout matrix and is marked unverified.
  compartmentData <- list(
    complex = list(
      analyte = "rilzabrutinib-BTK complex",
      units = "nmol/L",
      specimen = "whole blood",
      verified = FALSE
    )
  )

  covariateData <- list(
    CP_RILZABRUTINIB_NGML = list(
      description = "Time-varying TOTAL rilzabrutinib plasma concentration driving BTK binding",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "TOTAL (bound plus unbound) plasma concentration. The free",
        "concentration that drives Eq. 1 is computed inside model() as",
        "CP_RILZABRUTINIB_NGML * fup and converted to nmol/L with the",
        "molecular weight. In the source study this concentration came from",
        "a PK-Sim PBPK model that is not extracted. Published day-10",
        "predicted anchors (Yan 2026 Table 2): Cmax 192 ng/mL and AUC216-240",
        "466 ng*h/mL at 300 mg OD; Cmax 345 and AUC 1027 at 600 mg OD;",
        "Cmax 194 and AUC228-240 465 at 300 mg BID; Cmax 292 and AUC228-240",
        "696 at 450 mg BID. Set to 0 for drug-free records."
      ),
      source_name = "(none; plasma concentration predicted by the PK-Sim PBPK model)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 32L,
    n_studies = 1L,
    age_range = "23-27.5 years",
    weight_range = "72.8-91.4 kg",
    sex_female_pct = 0,
    race_ethnicity = "majority Caucasian",
    disease_state = "healthy volunteers (the paper applies the model to immune thrombocytopenia dosing)",
    dose_range = "rilzabrutinib 300 mg OD, 600 mg OD, 300 mg BID and 450 mg BID for 10 days (oral liquid formulation, fasted)",
    regions = NA_character_,
    notes = paste(
      "The BO data the model was optimised against came from the",
      "multiple-dose part (Part B) of the phase I study of Smith 2017",
      "(Br J Clin Pharmacol 83:2367; Yan 2026 ref 25), summarised in Yan",
      "2026 Supplementary Table S1: 32 healthy men, mean BMI 21.8-26.7",
      "kg/m2. BO was measured in whole blood (Table S1, Pharmacodynamic",
      "assessments). The PK-Sim virtual population (32 men) matched that",
      "design. The DDI and hepatic-impairment scenarios of the paper are",
      "PBPK-layer simulations and are not reproduced by this model, except",
      "that Eq. 2 changes the unbound fraction in hepatic impairment (see",
      "the fup parameter)."
    )
  )

  ini({
    mw_ril <- fixed(665.77)
    label("Molecular weight of rilzabrutinib (g/mol)")
    # Yan 2026 Table 1, Physicochemical and biochemical block, MW row
    # (665.77 g/mol). Converts the ng/mL covariate to the nmol/L scale of
    # kd and btk0.

    fup <- fixed(0.025)
    label("Fraction of rilzabrutinib unbound in plasma (unitless)")
    # Yan 2026 Table 1, PK parameters block, fup row (0.025, from the FDA
    # integrated review, ref 8). Healthy-subject value. For hepatic
    # impairment the paper scales plasma protein binding with Eq. 2,
    # PPSF = 1 / (fup + (1 - fup) * Albumin_f), so the impaired unbound
    # fraction is fup * PPSF; with the Supplementary Table S3 albumin
    # fractions 0.9 / 0.8 / 0.6 (mild / moderate / severe) that gives
    # 0.0277 / 0.0311 / 0.0410. Override fup in ini() to use them.

    kd <- fixed(1.3)
    label("Equilibrium dissociation constant of rilzabrutinib for BTK (nmol/L)")
    # Yan 2026 Table 1, BO parameters block, 'Ki against BTK (nM)' = 1.30
    # (from Owens 2022, ref 22). The paper calls it an inhibition constant,
    # but Eq. 1 uses koff/Ki as the association rate constant, which makes
    # it the equilibrium dissociation constant koff/kon; recorded as kd, as
    # in Zhang_2026_ribociclib_qsp.

    koff <- fixed(1.2e-5 * 3600)
    label("First-order dissociation rate constant of the RIL-BTK complex (1/h)")
    # Yan 2026 Table 1, BO parameters block, 'koff (x 10^-5 s^-1)' = 1.2,
    # 'Optimized using the value from Ref.24'. Converted to hours:
    # 1.2e-5 /s * 3600 s/h = 0.0432 /h (complex half-life 16 h). Methods:
    # koff 'was optimized to 1.2 x 10^-5 s^-1 from its original literature
    # value to improve the fit between observed and predicted BO'.

    btk0 <- fixed(100)
    label("Initial (total) BTK concentration (nmol/L)")
    # Yan 2026 Table 1, BO parameters block, 'Initial BTK (uM)' = 0.10
    # (from Xu 2022, ref 15); Methods: 'set to a default value of 0.1 uM'.
    # 0.1 uM = 100 nmol/L.
  })

  model({
    # Free rilzabrutinib driving BTK binding. Eq. 1 calls C_RIL 'the free
    # RIL concentration' without naming the matrix; free plasma is used.
    #   C[nmol/L] = C[ng/mL] / MW[g/mol] * 1000
    cu_nm <- CP_RILZABRUTINIB_NGML * fup / mw_ril * 1000

    # Yan 2026 Eq. 1:
    #   dN/dt = koff/Ki * BTK_free * C_RIL - koff * BTK_complex
    # with BTK_free = BTK_0 - BTK_complex (BTK is not resynthesised or
    # degraded in Eq. 1, so the total pool is constant).
    kon <- koff / kd
    btk_free <- btk0 - complex
    d/dt(complex) <- kon * btk_free * cu_nm - koff * complex

    # BTK occupancy: Yan 2026 Methods, 'the percentage of BTK targets bound
    # by RIL to form the drug-target complex'. occ is the fraction, bo the
    # percentage the paper reports.
    occ <- complex / btk0
    bo <- 100 * occ
  })
}
