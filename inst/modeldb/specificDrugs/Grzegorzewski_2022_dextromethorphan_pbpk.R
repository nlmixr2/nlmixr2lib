Grzegorzewski_2022_dextromethorphan_pbpk <- function() {
  description <- paste(
    "PBPK (whole-body, 76 ODEs, SBML/libroadrunner). Dextromethorphan",
    "(DXM), its CYP2D6 metabolite dextrorphan (DXO) and dextrorphan",
    "O-glucuronide (DXO-Glu) in healthy adults, built by Grzegorzewski,",
    "Brandhorst and Koenig (2022) from 36 curated clinical studies to study",
    "CYP2D6 metabolic phenotyping with the urinary cumulative metabolic",
    "ratio UCMR = DXM / (DXO + DXO-Glu). Oral DXM dissolves into the gut",
    "lumen; 55% of the dissolved dose is absorbed into gut tissue (where",
    "CYP3A4 N-demethylation removes part of it) and the remainder goes to",
    "feces. Gut, pancreas and spleen drain to the portal vein; liver,",
    "kidney, lung, forearm and a lumped rest-of-body compartment exchange",
    "with their plasma spaces by tissue/plasma partition coefficients, and",
    "plasma is sampled at the forearm (median cubital) vein. In the liver",
    "DXM is O-demethylated to DXO by CYP2D6 and CYP3A4 (irreversible",
    "Michaelis-Menten) and DXO is glucuronidated by UGT; all three",
    "species are excreted in urine by first-order renal processes. The",
    "CYP2D6 activity score (covariate CYP2D6, 0-4) scales the CYP2D6 Vmax",
    "proportionally and the Km by AS^-0.4; AS = 0 removes the CYP2D6",
    "pathway. Between-subject variability is a correlated log-normal",
    "distribution of hepatic CYP2D6 and CYP3A4 Km and Vmax derived from",
    "human liver microsome data; its dispersion was digitised from",
    "Figure 2 because the paper prints no numeric values. Another",
    "dextromethorphan model is available: modellib('TerHeine_2014_dextromethorphan')."
  )
  reference <- paste(
    "Grzegorzewski J, Brandhorst J, Koenig M. Physiologically based",
    "pharmacokinetic (PBPK) modeling of the role of CYP2D6 polymorphism",
    "for metabolic phenotyping with dextromethorphan. Front Pharmacol.",
    "2022;13:1029073. doi:10.3389/fphar.2022.1029073.",
    "Parameter values from Table 2; CYP2D6 activity-score means from",
    "Figure 2B; between-subject variability digitised from Figure 2A",
    "and 2B. The paper does not print its rate laws; they, and the",
    "physiological constants absent from Table 2, were taken from the",
    "model archive the paper names as the version used (Section 2.2):",
    "Grzegorzewski J, Koenig M. Physiologically based pharmacokinetic",
    "(PBPK) model of dextromethorphan v0.9.5. Zenodo.",
    "doi:10.5281/zenodo.7025683 (file models/dextromethorphan_body_flat.xml).",
    "See the vignette Errata for the dissolution-rate unit and the",
    "forearm-outflow reaction.",
    sep = " "
  )
  vignette <- "Grzegorzewski_2022_dextromethorphan_pbpk"
  units <- list(time = "min", dosing = "mg", concentration = "nmol/L")

  # Oral doses go into `depot` (mg DXM); intravenous doses into `depot_iv`
  # (mg DXM), which empties into venous plasma with a half-life equal to the
  # injection time (source SBML reaction iv_dxm).
  dosing <- c("depot", "depot_iv")

  # Every plasma and tissue state is a CONCENTRATION (mmol/L), following the
  # source SBML. `depot` / `depot_iv` hold mg; urine and feces hold mmol.
  paper_specific_compartments <- c(
    "feces",
    "hepatic_vein",
    "kidney_vein",
    "rest_plasma",
    "rest",
    "rest_vein",
    "forearm_plasma",
    "forearm",
    "forearm_vein",
    "hepatic_vein_dxor",
    "kidney_vein_dxor",
    "rest_plasma_dxor",
    "rest_dxor",
    "rest_vein_dxor",
    "forearm_plasma_dxor",
    "forearm_dxor",
    "forearm_vein_dxor",
    "hepatic_vein_dxorgluc",
    "kidney_vein_dxorgluc",
    "rest_plasma_dxorgluc",
    "rest_dxorgluc",
    "rest_vein_dxorgluc",
    "forearm_plasma_dxorgluc",
    "forearm_dxorgluc",
    "forearm_vein_dxorgluc"
  )

  compartmentData <- list(
    depot = list(analyte = "dextromethorphan", units = "mg", specimen = "administration site", verified = TRUE),
    depot_iv = list(analyte = "dextromethorphan", units = "mg", specimen = "administration site", verified = TRUE),
    gut_lumen = list(analyte = "dextromethorphan", units = "mmol/L", specimen = "administration site", verified = TRUE),
    feces = list(analyte = "dextromethorphan", units = "mmol", specimen = "faeces", verified = TRUE),
    gut_plasma = list(analyte = "dextromethorphan", units = "mmol/L", specimen = "plasma", verified = TRUE),
    gut = list(analyte = "dextromethorphan", units = "mmol/L", specimen = "tissue", verified = TRUE),
    pancreas_plasma = list(analyte = "dextromethorphan", units = "mmol/L", specimen = "plasma", verified = TRUE),
    pancreas = list(analyte = "dextromethorphan", units = "mmol/L", specimen = "tissue", verified = TRUE),
    spleen_plasma = list(analyte = "dextromethorphan", units = "mmol/L", specimen = "plasma", verified = TRUE),
    spleen = list(analyte = "dextromethorphan", units = "mmol/L", specimen = "tissue", verified = TRUE),
    portal = list(analyte = "dextromethorphan", units = "mmol/L", specimen = "plasma", verified = TRUE),
    liver_plasma = list(analyte = "dextromethorphan", units = "mmol/L", specimen = "plasma", verified = TRUE),
    liver = list(analyte = "dextromethorphan", units = "mmol/L", specimen = "tissue", verified = TRUE),
    hepatic_vein = list(analyte = "dextromethorphan", units = "mmol/L", specimen = "plasma", verified = TRUE),
    kidney_plasma = list(analyte = "dextromethorphan", units = "mmol/L", specimen = "plasma", verified = TRUE),
    kidney = list(analyte = "dextromethorphan", units = "mmol/L", specimen = "tissue", verified = TRUE),
    kidney_vein = list(analyte = "dextromethorphan", units = "mmol/L", specimen = "plasma", verified = TRUE),
    urine = list(analyte = "dextromethorphan", units = "mmol", specimen = "urine", verified = TRUE),
    lung_plasma = list(analyte = "dextromethorphan", units = "mmol/L", specimen = "plasma", verified = TRUE),
    lung = list(analyte = "dextromethorphan", units = "mmol/L", specimen = "tissue", verified = TRUE),
    rest_plasma = list(analyte = "dextromethorphan", units = "mmol/L", specimen = "plasma", verified = TRUE),
    rest = list(analyte = "dextromethorphan", units = "mmol/L", specimen = "tissue", verified = TRUE),
    rest_vein = list(analyte = "dextromethorphan", units = "mmol/L", specimen = "plasma", verified = TRUE),
    forearm_plasma = list(analyte = "dextromethorphan", units = "mmol/L", specimen = "plasma", verified = TRUE),
    forearm = list(analyte = "dextromethorphan", units = "mmol/L", specimen = "tissue", verified = TRUE),
    forearm_vein = list(analyte = "dextromethorphan", units = "mmol/L", specimen = "plasma", verified = TRUE),
    venous = list(analyte = "dextromethorphan", units = "mmol/L", specimen = "plasma", verified = TRUE),
    arterial = list(analyte = "dextromethorphan", units = "mmol/L", specimen = "plasma", verified = TRUE),
    gut_plasma_dxor = list(analyte = "dextrorphan", units = "mmol/L", specimen = "plasma", verified = TRUE),
    gut_dxor = list(analyte = "dextrorphan", units = "mmol/L", specimen = "tissue", verified = TRUE),
    pancreas_plasma_dxor = list(analyte = "dextrorphan", units = "mmol/L", specimen = "plasma", verified = TRUE),
    pancreas_dxor = list(analyte = "dextrorphan", units = "mmol/L", specimen = "tissue", verified = TRUE),
    spleen_plasma_dxor = list(analyte = "dextrorphan", units = "mmol/L", specimen = "plasma", verified = TRUE),
    spleen_dxor = list(analyte = "dextrorphan", units = "mmol/L", specimen = "tissue", verified = TRUE),
    portal_dxor = list(analyte = "dextrorphan", units = "mmol/L", specimen = "plasma", verified = TRUE),
    liver_plasma_dxor = list(analyte = "dextrorphan", units = "mmol/L", specimen = "plasma", verified = TRUE),
    liver_dxor = list(analyte = "dextrorphan", units = "mmol/L", specimen = "tissue", verified = TRUE),
    hepatic_vein_dxor = list(analyte = "dextrorphan", units = "mmol/L", specimen = "plasma", verified = TRUE),
    kidney_plasma_dxor = list(analyte = "dextrorphan", units = "mmol/L", specimen = "plasma", verified = TRUE),
    kidney_dxor = list(analyte = "dextrorphan", units = "mmol/L", specimen = "tissue", verified = TRUE),
    kidney_vein_dxor = list(analyte = "dextrorphan", units = "mmol/L", specimen = "plasma", verified = TRUE),
    urine_dxor = list(analyte = "dextrorphan", units = "mmol", specimen = "urine", verified = TRUE),
    lung_plasma_dxor = list(analyte = "dextrorphan", units = "mmol/L", specimen = "plasma", verified = TRUE),
    lung_dxor = list(analyte = "dextrorphan", units = "mmol/L", specimen = "tissue", verified = TRUE),
    rest_plasma_dxor = list(analyte = "dextrorphan", units = "mmol/L", specimen = "plasma", verified = TRUE),
    rest_dxor = list(analyte = "dextrorphan", units = "mmol/L", specimen = "tissue", verified = TRUE),
    rest_vein_dxor = list(analyte = "dextrorphan", units = "mmol/L", specimen = "plasma", verified = TRUE),
    forearm_plasma_dxor = list(analyte = "dextrorphan", units = "mmol/L", specimen = "plasma", verified = TRUE),
    forearm_dxor = list(analyte = "dextrorphan", units = "mmol/L", specimen = "tissue", verified = TRUE),
    forearm_vein_dxor = list(analyte = "dextrorphan", units = "mmol/L", specimen = "plasma", verified = TRUE),
    venous_dxor = list(analyte = "dextrorphan", units = "mmol/L", specimen = "plasma", verified = TRUE),
    arterial_dxor = list(analyte = "dextrorphan", units = "mmol/L", specimen = "plasma", verified = TRUE),
    gut_plasma_dxorgluc = list(
      analyte = "dextrorphan O-glucuronide",
      units = "mmol/L",
      specimen = "plasma",
      verified = TRUE
    ),
    gut_dxorgluc = list(analyte = "dextrorphan O-glucuronide", units = "mmol/L", specimen = "tissue", verified = TRUE),
    pancreas_plasma_dxorgluc = list(
      analyte = "dextrorphan O-glucuronide",
      units = "mmol/L",
      specimen = "plasma",
      verified = TRUE
    ),
    pancreas_dxorgluc = list(
      analyte = "dextrorphan O-glucuronide",
      units = "mmol/L",
      specimen = "tissue",
      verified = TRUE
    ),
    spleen_plasma_dxorgluc = list(
      analyte = "dextrorphan O-glucuronide",
      units = "mmol/L",
      specimen = "plasma",
      verified = TRUE
    ),
    spleen_dxorgluc = list(
      analyte = "dextrorphan O-glucuronide",
      units = "mmol/L",
      specimen = "tissue",
      verified = TRUE
    ),
    portal_dxorgluc = list(
      analyte = "dextrorphan O-glucuronide",
      units = "mmol/L",
      specimen = "plasma",
      verified = TRUE
    ),
    liver_plasma_dxorgluc = list(
      analyte = "dextrorphan O-glucuronide",
      units = "mmol/L",
      specimen = "plasma",
      verified = TRUE
    ),
    liver_dxorgluc = list(
      analyte = "dextrorphan O-glucuronide",
      units = "mmol/L",
      specimen = "tissue",
      verified = TRUE
    ),
    hepatic_vein_dxorgluc = list(
      analyte = "dextrorphan O-glucuronide",
      units = "mmol/L",
      specimen = "plasma",
      verified = TRUE
    ),
    kidney_plasma_dxorgluc = list(
      analyte = "dextrorphan O-glucuronide",
      units = "mmol/L",
      specimen = "plasma",
      verified = TRUE
    ),
    kidney_dxorgluc = list(
      analyte = "dextrorphan O-glucuronide",
      units = "mmol/L",
      specimen = "tissue",
      verified = TRUE
    ),
    kidney_vein_dxorgluc = list(
      analyte = "dextrorphan O-glucuronide",
      units = "mmol/L",
      specimen = "plasma",
      verified = TRUE
    ),
    urine_dxorgluc = list(analyte = "dextrorphan O-glucuronide", units = "mmol", specimen = "urine", verified = TRUE),
    lung_plasma_dxorgluc = list(
      analyte = "dextrorphan O-glucuronide",
      units = "mmol/L",
      specimen = "plasma",
      verified = TRUE
    ),
    lung_dxorgluc = list(analyte = "dextrorphan O-glucuronide", units = "mmol/L", specimen = "tissue", verified = TRUE),
    rest_plasma_dxorgluc = list(
      analyte = "dextrorphan O-glucuronide",
      units = "mmol/L",
      specimen = "plasma",
      verified = TRUE
    ),
    rest_dxorgluc = list(analyte = "dextrorphan O-glucuronide", units = "mmol/L", specimen = "tissue", verified = TRUE),
    rest_vein_dxorgluc = list(
      analyte = "dextrorphan O-glucuronide",
      units = "mmol/L",
      specimen = "plasma",
      verified = TRUE
    ),
    forearm_plasma_dxorgluc = list(
      analyte = "dextrorphan O-glucuronide",
      units = "mmol/L",
      specimen = "plasma",
      verified = TRUE
    ),
    forearm_dxorgluc = list(
      analyte = "dextrorphan O-glucuronide",
      units = "mmol/L",
      specimen = "tissue",
      verified = TRUE
    ),
    forearm_vein_dxorgluc = list(
      analyte = "dextrorphan O-glucuronide",
      units = "mmol/L",
      specimen = "plasma",
      verified = TRUE
    ),
    venous_dxorgluc = list(
      analyte = "dextrorphan O-glucuronide",
      units = "mmol/L",
      specimen = "plasma",
      verified = TRUE
    ),
    arterial_dxorgluc = list(
      analyte = "dextrorphan O-glucuronide",
      units = "mmol/L",
      specimen = "plasma",
      verified = TRUE
    )
  )
  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Scales every absolute organ volume (fractional volumes in L/kg",
        "times body weight) and, through cardiac output per body weight",
        "(COBW = 1.548 mL/s/kg), every blood flow. The reference",
        "individual is 75 kg (Table 2, BW, ICRP male). The forearm tissue",
        "and forearm-vein volumes are fixed at 1 L in the source archive",
        "and do not scale with body weight."
      ),
      source_name = "BW"
    ),
    CYP2D6 = list(
      description = "CYP2D6 activity score (AS), sum of the two allele activity values",
      units = "(dimensionless)",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Activity score from the CYP2D6 diplotype (PharmGKB allele",
        "activities, Supplementary Table S1: *1, *2 = 1; *10 = 0.25;",
        "*17, *29, *41 = 0.5; *3, *4, *5, *6 = 0; gene duplications add",
        "the duplicated activity). Range 0-4 in the paper (Table 2",
        "LI__cyp2d6_ac '0.0 - 3.0'; Supplementary Table S2 also lists AS",
        "4). The CYP2D6 Vmax is proportional to AS and the Km is",
        "multiplied by AS^lambda_1 with lambda_1 = -0.4 (Table 2),",
        "reproducing the Figure 2B means (Km 7.9 uM at AS 1, 6.0 uM at",
        "AS 2). AS = 0 means no CYP2D6 activity, leaving only the minor",
        "CYP3A4 O-demethylation. There is no reference value; the source",
        "archive defaults to AS = 2 (*1/*1). Source name LI__cyp2d6_ac",
        "(the archive also accepts two allele indices, which map to the",
        "same activity sum)."
      ),
      source_name = "LI__cyp2d6_ac"
    )
  )

  population <- list(
    species = "human",
    n_subjects = NA_integer_,
    n_studies = 36L,
    age_range = "adults (age >= 18 years, Section 2.1 eligibility)",
    weight_range = "75 kg reference individual (Table 2); per-study weights not reported",
    sex_female_pct = NA_real_,
    race_ethnicity = paste(
      "Curated studies span European, East Asian, Mexican Mestizo, Cuban",
      "and Trinidadian cohorts; population simulations use the",
      "biogeographical AS frequencies of Supplementary Table S2."
    ),
    disease_state = paste(
      "Healthy volunteers; studies with drug-drug interactions known to",
      "alter DXM pharmacokinetics were excluded (Section 2.1)."
    ),
    dose_range = paste(
      "Oral 2 mg to 3 mg/kg DXM or DXM hydrobromide (solution, syrup,",
      "capsule, tablet, sustained release) and one intravenous study",
      "(Duedahl 2005, 0.5 mg/kg); most studies gave 30 mg (Table 1)."
    ),
    regions = "Multinational (Table 1)",
    notes = paste(
      "PBPK model calibrated against time-courses from 36 studies curated",
      "into PK-DB (Table 1) out of 404 screened records. A subset of",
      "parameters (flagged F in Table 2) was fitted to the data of Figures",
      "4-9; the rest are physiological reference values. Variability in",
      "CYP2D6 and CYP3A4 Km and Vmax was fitted to human liver microsome",
      "data (Yang 2012; Storelli 2019a), not to the clinical data. The",
      "CYP2D6 activity-score frequencies of the curated UCMR data",
      "(Figure 2B, P(AS)) were 0.06, 0.02, 0.07, 0.29, 0.13, 0.06, 0.33",
      "and 0.02 for AS 0, 0.25, 0.5, 1, 1.25, 1.5, 2 and 3."
    )
  )

  ini({
    # ============ Reference physiology (Table 2; SBML archive) ============
    # The whole-body circulation ODEs are not printed in the paper; they and
    # the physiological constants Table 2 omits were taken from the cited
    # model archive (Zenodo 7025683, dextromethorphan_body_flat.xml).
    cobw <- fixed(1.548); label("Cardiac output per body weight (mL/s/kg)")  # Table 2 COBW
    hct <- fixed(0.51); label("Hematocrit (fraction)")  # Table 2 HCT, upper range male
    fblood <- fixed(0.02); label("Blood fraction of organ volume (unitless)")  # SBML Fblood

    # Fractional tissue volumes (L/kg); forearm is fixed at 1 L in the
    # archive and so does not scale with body weight.
    fvgu <- fixed(0.0171); label("Gut fractional tissue volume (L/kg)")  # Table 2 FVgu
    fvki <- fixed(0.0044); label("Kidney fractional tissue volume (L/kg)")  # Table 2 FVki
    fvli <- fixed(0.021); label("Liver fractional tissue volume (L/kg)")  # Table 2 FVli
    fvlu <- fixed(0.0076); label("Lung fractional tissue volume (L/kg)")  # Table 2 FVlu
    fvsp <- fixed(0.0026); label("Spleen fractional tissue volume (L/kg)")  # Table 2 FVsp
    fvpa <- fixed(0.01); label("Pancreas fractional tissue volume (L/kg)")  # Table 2 FVpa
    fvfo <- fixed(0.00482857142857143); label("Forearm fractional volume used only in the rest-of-body balance (L/kg)")  # Table 2 FVfo; archive FVfo
    fvve <- fixed(0.0514); label("Venous blood fractional volume (L/kg)")  # Table 2 FVve
    fvar <- fixed(0.0257); label("Arterial blood fractional volume (L/kg)")  # Table 2 FVar
    fvpo <- fixed(0.001); label("Portal plasma fractional volume (L/kg)")  # Table 2 FVpo
    fvhv <- fixed(0.001); label("Hepatic venous plasma fractional volume (L/kg)")  # SBML FVhv
    fvrev <- fixed(0.001); label("Rest venous plasma fractional volume (L/kg)")  # SBML FVrev
    fvkiv <- fixed(0.001); label("Kidney venous plasma fractional volume (L/kg)")  # SBML FVkiv
    fvfov <- fixed(0.001); label("Forearm venous plasma fractional volume (L/kg)")  # SBML FVfov

    # Fractional tissue blood flows (fraction of cardiac output)
    fqgu <- fixed(0.146); label("Gut fractional blood flow (unitless)")  # Table 2 FQgu
    fqki <- fixed(0.19); label("Kidney fractional blood flow (unitless)")  # Table 2 FQki
    fqh <- fixed(0.215); label("Hepatic venous-side fractional blood flow (unitless)")  # Table 2 FQh
    fqlu <- fixed(1); label("Lung fractional blood flow (unitless)")  # Table 2 FQlu
    fqsp <- fixed(0.017); label("Spleen fractional blood flow (unitless)")  # Table 2 FQsp
    fqpa <- fixed(0.017); label("Pancreas fractional blood flow (unitless)")  # Table 2 FQpa
    fqfo <- fixed(0.0146153846153846); label("Forearm fractional blood flow (unitless)")  # Table 2 FQfo
    f_shunting_forearm <- fixed(0.279517545017716); label("Fraction of forearm arterial flow shunted to the forearm vein (unitless)")  # Table 2 f_shunting_forearm (fitted)

    # ============ Molecular weight (Table 2) ============
    # Only the parent MW is needed: it converts the mg dose to mmol at
    # dissolution and intravenous release. The states and observations are
    # concentrations (mM / nmol/L), so the metabolite molecular weights
    # Mr_dxo and Mr_dxo_glu (Table 2, for mg-unit outputs) are not required.
    mr_dxm <- fixed(271.404); label("Molecular weight of dextromethorphan (g/mol)")  # SBML Mr_dxm

    # ============ Distribution (Table 2, fitted) ============
    kp_dxm <- fixed(8.73461301162365); label("Tissue/plasma partition coefficient, dextromethorphan (unitless)")  # Table 2 Kp_dxm = 8.7346 (fitted)
    kp_fo_dxm <- fixed(10); label("Tissue/plasma partition coefficient, dextromethorphan in forearm (unitless)")  # Table 2 Kp_fo_dxm = 10 (fitted)
    kp_dxo <- fixed(4); label("Tissue/plasma partition coefficient, dextrorphan (unitless)")  # Table 2 Kp_dxo = 4 (fitted)
    kp_dxoglu <- fixed(0.08); label("Tissue/plasma partition coefficient, dextrorphan O-glucuronide (unitless)")  # Table 2 Kp_dxo_glu = 0.08 (fitted)
    ftissue_dxm <- fixed(1000); label("Tissue distribution rate, dextromethorphan (L/min)")  # Table 2 ftissue_dxm (fitted)
    ftissue_dxo <- fixed(100); label("Tissue distribution rate, dextrorphan (L/min)")  # Table 2 ftissue_dxo (fitted)
    ftissue_dxoglu <- fixed(3); label("Tissue distribution rate, dextrorphan O-glucuronide (L/min)")  # Table 2 ftissue_dxo_glu (fitted)

    # ============ Absorption (Table 2) ============
    # The archive stores the dissolution rate as 0.0217 /hr, but the Table 2
    # header reads 1/hr while a 0.0217/hr dissolution gives a >24 h Tmax that
    # contradicts the 1-2 h Tmax of every oral panel in Figures 3-4. Treating
    # the rate as 0.0217/min reproduces those figures. See the vignette Errata.
    ka_dis_dxm <- fixed(1.30434782608696); label("Dissolution + stomach passage rate, dextromethorphan (1/hr)")  # Table 2 Ka_dis_dxm 0.0217 reinterpreted as 1/min = 1.304/hr
    ka_abs_dxm <- fixed(3.42854395495006); label("Intestinal absorption rate, dextromethorphan (1/hr)")  # Table 2 GU__Ka_abs_dxm (fitted)
    f_dxm_abs <- fixed(0.55); label("Fraction of dissolved dextromethorphan absorbed (unitless)")  # Table 2 GU__F_dxm = 0.55 (Schadel 1995)

    # ============ Intestinal first-pass CYP3A4 (Table 2) ============
    gu_cyp3a4_vmax <- fixed(0.0002); label("Vmax, gut CYP3A4 N-demethylation of dextromethorphan (mmol/min/L)")  # Table 2 GU__DXMCYP3A4_Vmax (fitted)
    gu_cyp3a4_km <- fixed(0.7); label("Km, gut CYP3A4 N-demethylation of dextromethorphan (mmol/L)")  # Table 2 GU__DXMCYP3A4_Km (Kerry 1994; Yu 2001)

    # ============ Hepatic metabolism (Table 2) ============
    # The CYP2D6 and CYP3A4 Km and Vmax are log-scale typical values (fixed to
    # the Table 2 point estimates) carrying between-subject variability; see
    # the variability block below.
    lvmax2d6 <- fixed(log(0.003)); label("Vmax, hepatic CYP2D6 O-demethylation of dextromethorphan (log mmol/min/L)")  # Table 2 LI__DXMCYP2D6_Vmax (fitted)
    lkm2d6 <- fixed(log(0.0079)); label("Km, hepatic CYP2D6 O-demethylation of dextromethorphan (log mmol/L)")  # Table 2 LI__DXMCYP2D6_Km (Storelli 2019a; Yang 2012)
    cyp2d6_lambda1 <- fixed(-0.4); label("Activity-score scaling exponent for CYP2D6 Km (unitless)")  # Table 2 LI__lambda_1 (fitted)
    lvmax3a4 <- fixed(log(0.0004)); label("Vmax, hepatic CYP3A4 O-demethylation of dextromethorphan (log mmol/min/L)")  # Table 2 LI__DXMCYP3A4_Vmax (fitted)
    lkm3a4 <- fixed(log(0.157)); label("Km, hepatic CYP3A4 O-demethylation of dextromethorphan (log mmol/L)")  # Table 2 LI__DXMCYP3A4_Km (Yu 2001)
    li_ugt_vmax <- fixed(0.895285026249854); label("Vmax, hepatic UGT glucuronidation of dextrorphan (mmol/min/L)")  # Table 2 LI__DXOUGT_Vmax (fitted)
    li_ugt_km <- fixed(0.69); label("Km, hepatic UGT glucuronidation of dextrorphan (mmol/L)")  # Table 2 LI__DXOUGT_Km (Lutz 2012)

    # ============ Renal excretion (Table 2, fitted) ============
    kex_dxm <- fixed(0.017); label("Urinary excretion rate, dextromethorphan (1/min)")  # Table 2 KI__DXMEX_k (fitted)
    kex_dxo <- fixed(0.3); label("Urinary excretion rate, dextrorphan (1/min)")  # Table 2 KI__DXOEX_k (fitted)
    kex_dxoglu <- fixed(10); label("Urinary excretion rate, dextrorphan O-glucuronide (1/min)")  # Table 2 KI__DXOGLUEX_k (fitted)

    # ============ Intravenous input ============
    ti_dxm <- fixed(10); label("Nominal intravenous injection time (s)")  # SBML ti_dxm; first-order release with half-life ti_dxm

    # ============ Between-subject variability ============
    # Correlated log-normal distributions of hepatic CYP2D6 and CYP3A4 Km and
    # Vmax, fitted to human liver microsome data (Section 2.3). The paper
    # reports no numeric dispersion, so the variances and correlations below
    # were digitised from Figure 2 (CYP3A4 from the normalised panel 2A;
    # CYP2D6 within-activity-score from the panels of 2B). CYP2D6 and CYP3A4
    # are independent. See the vignette Errata.
    etalvmax2d6 + etalkm2d6 ~ c(
      0.85,
      -0.092, 0.44
    )
    etalvmax3a4 + etalkm3a4 ~ c(
      0.69,
      0.18, 0.45
    )

    # ============ Residual error ============
    propSd <- fixed(0); label("Proportional residual error on dextromethorphan (fraction; ZERO - not reported by the paper)")  # Section 2.3: variability is modelled through the enzyme Km/Vmax distributions, not a residual-error term
  })

  model({
    # ============ Scenario / covariate ============
    # CYP2D6 activity score (0-4). Vmax scales with AS; Km scales with
    # AS^lambda_1; AS = 0 removes the CYP2D6 pathway (Section 2.3).
    cyp2d6_activity <- CYP2D6

    # ============ Organ volumes (L); forearm is fixed at 1 L ============
    fvre <- 1 - (fvgu + fvki + fvli + fvlu + fvsp + fvpa + fvve + fvar + fvfo)

    vgu <- WT * fvgu
    vki <- WT * fvki
    vli <- WT * fvli
    vlu <- WT * fvlu
    vsp <- WT * fvsp
    vpa <- WT * fvpa
    vre <- WT * fvre

    # Blood pools corrected for the blood held within the organ volumes.
    vve <- WT * fvve - (fvve / (fvar + fvve)) * WT * fblood * (1 - fvve - fvar)
    var_ <- WT * fvar - (fvar / (fvar + fvve)) * WT * fblood * (1 - fvve - fvar)
    vpo <- (1 - hct) * (WT * fvpo - (fvpo / (fvar + fvve + fvpo + fvhv + fvkiv + fvrev + fvfov)) * WT * fblood * (1 - (fvar + fvve + fvpo + fvhv + fvkiv + fvrev + fvfov)))
    vhv <- (1 - hct) * (WT * fvhv - (fvhv / (fvar + fvve + fvpo + fvhv + fvkiv + fvrev + fvfov)) * WT * fblood * (1 - (fvar + fvve + fvpo + fvhv + fvkiv + fvrev + fvfov)))
    vrev <- WT * fvrev - (fvrev / (fvar + fvve)) * WT * fblood * (1 - fvve - fvar)
    vkiv <- WT * fvkiv - (fvkiv / (fvar + fvve)) * WT * fblood * (1 - fvve - fvar)
    vfov <- 1

    vgu_tissue <- vgu * (1 - fblood)
    vki_tissue <- vki * (1 - fblood)
    vli_tissue <- vli * (1 - fblood)
    vlu_tissue <- vlu * (1 - fblood)
    vpa_tissue <- vpa * (1 - fblood)
    vsp_tissue <- vsp * (1 - fblood)
    vre_tissue <- vre * (1 - fblood)
    vfo_tissue <- 1 * (1 - fblood)

    vgu_plasma <- vgu * fblood * (1 - hct)
    vki_plasma <- vki * fblood * (1 - hct)
    vli_plasma <- vli * fblood * (1 - hct)
    vlu_plasma <- vlu * fblood * (1 - hct)
    vpa_plasma <- vpa * fblood * (1 - hct)
    vsp_plasma <- vsp * fblood * (1 - hct)
    vre_plasma <- vre * fblood * (1 - hct)
    vfo_plasma <- 1 * fblood * (1 - hct)

    # ============ Blood flows (L/min) ============
    co <- WT * cobw
    qc <- (co / 1000) * 60
    fqre <- 1 - (fqki + fqh + fqfo)

    qgu <- qc * fqgu
    qki <- qc * fqki
    qh <- qc * fqh
    qlu <- qc * fqlu
    qsp <- qc * fqsp
    qpa <- qc * fqpa
    qfo <- qc * fqfo
    qre <- qc * fqre
    qpo <- qsp + qpa + qgu
    qha <- qh - qgu - qsp - qpa

    # ============ Absorption (oral), SBML intestine model ============
    vgulumen <- 1
    dissolution <- (ka_dis_dxm / 60) * depot / mr_dxm
    absorption <- (vgu_tissue * ka_abs_dxm / 60) * gut_lumen

    # ============ Intravenous input ============
    ki_dxm <- (0.693 / ti_dxm) * 60
    iv_dxm <- ki_dxm * depot_iv / mr_dxm

    # ============ Metabolism, SBML liver + intestine models ============
    # Hepatic CYP2D6/CYP3A4 Km and Vmax carry the correlated log-normal IIV.
    vmax2d6 <- exp(lvmax2d6 + etalvmax2d6)
    km2d6 <- exp(lkm2d6 + etalkm2d6)
    vmax3a4 <- exp(lvmax3a4 + etalvmax3a4)
    km3a4 <- exp(lkm3a4 + etalkm3a4)

    met_li_cyp2d6 <- vmax2d6 * vli_tissue * cyp2d6_activity * liver /
      (cyp2d6_activity^cyp2d6_lambda1 * km2d6 + liver)
    met_li_cyp3a4 <- vmax3a4 * vli_tissue * liver / (km3a4 + liver)
    met_li_ugt <- li_ugt_vmax * vli_tissue * liver_dxor / (li_ugt_km + liver_dxor)
    met_gu_cyp3a4 <- gu_cyp3a4_vmax * vgu_tissue * gut / (gut + gu_cyp3a4_km)

    # ============ Dosing compartments ============
    # depot / depot_iv hold mg; dissolution and the iv release return mmol.
    d/dt(depot) <- -dissolution * mr_dxm
    d/dt(depot_iv) <- -iv_dxm * mr_dxm

    # Intestinal lumen and feces (SBML intestine model). The gut-lumen
    # compartment has a fixed 1 L volume in the archive, so concentration
    # equals amount; both the absorbed and the excreted fractions leave it.
    d/dt(gut_lumen) <- (dissolution - absorption) / vgulumen
    d/dt(feces) <- (1 - f_dxm_abs) * absorption

    # ----- dextromethorphan: tissue uptake (SBML transport_<organ>_dxm) and forearm mixing -----
    tr_gu <- ftissue_dxm * (gut_plasma * kp_dxm - gut)
    tr_pa <- ftissue_dxm * (pancreas_plasma * kp_dxm - pancreas)
    tr_sp <- ftissue_dxm * (spleen_plasma * kp_dxm - spleen)
    tr_li <- ftissue_dxm * (liver_plasma * kp_dxm - liver)
    tr_ki <- ftissue_dxm * (kidney_plasma * kp_dxm - kidney)
    tr_lu <- ftissue_dxm * (lung_plasma * kp_dxm - lung)
    tr_re <- ftissue_dxm * (rest_plasma * kp_dxm - rest)
    tr_fo <- ftissue_dxm * (forearm_plasma * kp_fo_dxm - forearm)
    # SBML Flow_fo_fov_dxm: forearm plasma AND arterial are both stoich-1
    # reactants, so each loses the full summed rate (see vignette Errata).
    flow_fo_fov <- qfo * (1 - f_shunting_forearm) * forearm_plasma + qfo * f_shunting_forearm * arterial

    # ----- dextrorphan: tissue uptake (SBML transport_<organ>_dxo) and forearm mixing -----
    tr_gu_dxor <- ftissue_dxo * (gut_plasma_dxor * kp_dxo - gut_dxor)
    tr_pa_dxor <- ftissue_dxo * (pancreas_plasma_dxor * kp_dxo - pancreas_dxor)
    tr_sp_dxor <- ftissue_dxo * (spleen_plasma_dxor * kp_dxo - spleen_dxor)
    tr_li_dxor <- ftissue_dxo * (liver_plasma_dxor * kp_dxo - liver_dxor)
    tr_ki_dxor <- ftissue_dxo * (kidney_plasma_dxor * kp_dxo - kidney_dxor)
    tr_lu_dxor <- ftissue_dxo * (lung_plasma_dxor * kp_dxo - lung_dxor)
    tr_re_dxor <- ftissue_dxo * (rest_plasma_dxor * kp_dxo - rest_dxor)
    tr_fo_dxor <- ftissue_dxo * (forearm_plasma_dxor * kp_dxo - forearm_dxor)
    # SBML Flow_fo_fov_dxo: forearm plasma AND arterial are both stoich-1
    # reactants, so each loses the full summed rate (see vignette Errata).
    flow_fo_fov_dxor <- qfo * (1 - f_shunting_forearm) * forearm_plasma_dxor + qfo * f_shunting_forearm * arterial_dxor

    # ----- dextrorphan O-glucuronide: tissue uptake (SBML transport_<organ>_dxoglu) and forearm mixing -----
    tr_gu_dxorgluc <- ftissue_dxoglu * (gut_plasma_dxorgluc * kp_dxoglu - gut_dxorgluc)
    tr_pa_dxorgluc <- ftissue_dxoglu * (pancreas_plasma_dxorgluc * kp_dxoglu - pancreas_dxorgluc)
    tr_sp_dxorgluc <- ftissue_dxoglu * (spleen_plasma_dxorgluc * kp_dxoglu - spleen_dxorgluc)
    tr_li_dxorgluc <- ftissue_dxoglu * (liver_plasma_dxorgluc * kp_dxoglu - liver_dxorgluc)
    tr_ki_dxorgluc <- ftissue_dxoglu * (kidney_plasma_dxorgluc * kp_dxoglu - kidney_dxorgluc)
    tr_lu_dxorgluc <- ftissue_dxoglu * (lung_plasma_dxorgluc * kp_dxoglu - lung_dxorgluc)
    tr_re_dxorgluc <- ftissue_dxoglu * (rest_plasma_dxorgluc * kp_dxoglu - rest_dxorgluc)
    tr_fo_dxorgluc <- ftissue_dxoglu * (forearm_plasma_dxorgluc * kp_dxoglu - forearm_dxorgluc)
    # SBML Flow_fo_fov_dxoglu: forearm plasma AND arterial are both stoich-1
    # reactants, so each loses the full summed rate (see vignette Errata).
    flow_fo_fov_dxorgluc <- qfo * (1 - f_shunting_forearm) * forearm_plasma_dxorgluc + qfo * f_shunting_forearm * arterial_dxorgluc

    # ----- ODEs: dextromethorphan -----
    d/dt(gut_plasma) <- (qgu * arterial - qgu * gut_plasma - tr_gu) / vgu_plasma
    d/dt(gut) <- (tr_gu + f_dxm_abs * absorption - met_gu_cyp3a4) / vgu_tissue
    d/dt(pancreas_plasma) <- (qpa * arterial - qpa * pancreas_plasma - tr_pa) / vpa_plasma
    d/dt(pancreas) <- tr_pa / vpa_tissue
    d/dt(spleen_plasma) <- (qsp * arterial - qsp * spleen_plasma - tr_sp) / vsp_plasma
    d/dt(spleen) <- tr_sp / vsp_tissue
    d/dt(portal) <- (qgu * gut_plasma + qpa * pancreas_plasma + qsp * spleen_plasma - qpo * portal) / vpo
    d/dt(liver_plasma) <- (qha * arterial + qpo * portal - qh * liver_plasma - tr_li) / vli_plasma
    d/dt(liver) <- (tr_li - met_li_cyp2d6 - met_li_cyp3a4) / vli_tissue
    d/dt(hepatic_vein) <- (qh * liver_plasma - qh * hepatic_vein) / vhv
    d/dt(kidney_plasma) <- (qki * arterial - qki * kidney_plasma - tr_ki) / vki_plasma
    d/dt(kidney) <- (tr_ki - kex_dxm * vki_tissue * kidney) / vki_tissue
    d/dt(kidney_vein) <- (qki * kidney_plasma - qki * kidney_vein) / vkiv
    d/dt(urine) <- kex_dxm * vki_tissue * kidney
    d/dt(lung_plasma) <- (qlu * venous - qlu * lung_plasma - tr_lu) / vlu_plasma
    d/dt(lung) <- tr_lu / vlu_tissue
    d/dt(rest_plasma) <- (qre * arterial - qre * rest_plasma - tr_re) / vre_plasma
    d/dt(rest) <- tr_re / vre_tissue
    d/dt(rest_vein) <- (qre * rest_plasma - qre * rest_vein) / vrev
    d/dt(forearm_plasma) <- (qfo * (1 - f_shunting_forearm) * arterial - flow_fo_fov - tr_fo) / vfo_plasma
    d/dt(forearm) <- tr_fo / vfo_tissue
    d/dt(forearm_vein) <- (flow_fo_fov - qfo * forearm_vein) / vfov
    d/dt(venous) <- (iv_dxm + qki * kidney_vein + qh * hepatic_vein + qre * rest_vein + qfo * forearm_vein - qlu * venous) / vve
    d/dt(arterial) <- (qlu * lung_plasma - (qgu + qki + qha + qpa + qsp + qre + qfo * (1 - f_shunting_forearm)) * arterial - flow_fo_fov) / var_

    # ----- ODEs: dextrorphan -----
    d/dt(gut_plasma_dxor) <- (qgu * arterial_dxor - qgu * gut_plasma_dxor - tr_gu_dxor) / vgu_plasma
    d/dt(gut_dxor) <- (tr_gu_dxor) / vgu_tissue
    d/dt(pancreas_plasma_dxor) <- (qpa * arterial_dxor - qpa * pancreas_plasma_dxor - tr_pa_dxor) / vpa_plasma
    d/dt(pancreas_dxor) <- tr_pa_dxor / vpa_tissue
    d/dt(spleen_plasma_dxor) <- (qsp * arterial_dxor - qsp * spleen_plasma_dxor - tr_sp_dxor) / vsp_plasma
    d/dt(spleen_dxor) <- tr_sp_dxor / vsp_tissue
    d/dt(portal_dxor) <- (qgu * gut_plasma_dxor + qpa * pancreas_plasma_dxor + qsp * spleen_plasma_dxor - qpo * portal_dxor) / vpo
    d/dt(liver_plasma_dxor) <- (qha * arterial_dxor + qpo * portal_dxor - qh * liver_plasma_dxor - tr_li_dxor) / vli_plasma
    d/dt(liver_dxor) <- (tr_li_dxor + met_li_cyp2d6 + met_li_cyp3a4 - met_li_ugt) / vli_tissue
    d/dt(hepatic_vein_dxor) <- (qh * liver_plasma_dxor - qh * hepatic_vein_dxor) / vhv
    d/dt(kidney_plasma_dxor) <- (qki * arterial_dxor - qki * kidney_plasma_dxor - tr_ki_dxor) / vki_plasma
    d/dt(kidney_dxor) <- (tr_ki_dxor - kex_dxo * vki_tissue * kidney_dxor) / vki_tissue
    d/dt(kidney_vein_dxor) <- (qki * kidney_plasma_dxor - qki * kidney_vein_dxor) / vkiv
    d/dt(urine_dxor) <- kex_dxo * vki_tissue * kidney_dxor
    d/dt(lung_plasma_dxor) <- (qlu * venous_dxor - qlu * lung_plasma_dxor - tr_lu_dxor) / vlu_plasma
    d/dt(lung_dxor) <- tr_lu_dxor / vlu_tissue
    d/dt(rest_plasma_dxor) <- (qre * arterial_dxor - qre * rest_plasma_dxor - tr_re_dxor) / vre_plasma
    d/dt(rest_dxor) <- tr_re_dxor / vre_tissue
    d/dt(rest_vein_dxor) <- (qre * rest_plasma_dxor - qre * rest_vein_dxor) / vrev
    d/dt(forearm_plasma_dxor) <- (qfo * (1 - f_shunting_forearm) * arterial_dxor - flow_fo_fov_dxor - tr_fo_dxor) / vfo_plasma
    d/dt(forearm_dxor) <- tr_fo_dxor / vfo_tissue
    d/dt(forearm_vein_dxor) <- (flow_fo_fov_dxor - qfo * forearm_vein_dxor) / vfov
    d/dt(venous_dxor) <- (qki * kidney_vein_dxor + qh * hepatic_vein_dxor + qre * rest_vein_dxor + qfo * forearm_vein_dxor - qlu * venous_dxor) / vve
    d/dt(arterial_dxor) <- (qlu * lung_plasma_dxor - (qgu + qki + qha + qpa + qsp + qre + qfo * (1 - f_shunting_forearm)) * arterial_dxor - flow_fo_fov_dxor) / var_

    # ----- ODEs: dextrorphan O-glucuronide -----
    d/dt(gut_plasma_dxorgluc) <- (qgu * arterial_dxorgluc - qgu * gut_plasma_dxorgluc - tr_gu_dxorgluc) / vgu_plasma
    d/dt(gut_dxorgluc) <- (tr_gu_dxorgluc) / vgu_tissue
    d/dt(pancreas_plasma_dxorgluc) <- (qpa * arterial_dxorgluc - qpa * pancreas_plasma_dxorgluc - tr_pa_dxorgluc) / vpa_plasma
    d/dt(pancreas_dxorgluc) <- tr_pa_dxorgluc / vpa_tissue
    d/dt(spleen_plasma_dxorgluc) <- (qsp * arterial_dxorgluc - qsp * spleen_plasma_dxorgluc - tr_sp_dxorgluc) / vsp_plasma
    d/dt(spleen_dxorgluc) <- tr_sp_dxorgluc / vsp_tissue
    d/dt(portal_dxorgluc) <- (qgu * gut_plasma_dxorgluc + qpa * pancreas_plasma_dxorgluc + qsp * spleen_plasma_dxorgluc - qpo * portal_dxorgluc) / vpo
    d/dt(liver_plasma_dxorgluc) <- (qha * arterial_dxorgluc + qpo * portal_dxorgluc - qh * liver_plasma_dxorgluc - tr_li_dxorgluc) / vli_plasma
    d/dt(liver_dxorgluc) <- (tr_li_dxorgluc + met_li_ugt) / vli_tissue
    d/dt(hepatic_vein_dxorgluc) <- (qh * liver_plasma_dxorgluc - qh * hepatic_vein_dxorgluc) / vhv
    d/dt(kidney_plasma_dxorgluc) <- (qki * arterial_dxorgluc - qki * kidney_plasma_dxorgluc - tr_ki_dxorgluc) / vki_plasma
    d/dt(kidney_dxorgluc) <- (tr_ki_dxorgluc - kex_dxoglu * vki_tissue * kidney_dxorgluc) / vki_tissue
    d/dt(kidney_vein_dxorgluc) <- (qki * kidney_plasma_dxorgluc - qki * kidney_vein_dxorgluc) / vkiv
    d/dt(urine_dxorgluc) <- kex_dxoglu * vki_tissue * kidney_dxorgluc
    d/dt(lung_plasma_dxorgluc) <- (qlu * venous_dxorgluc - qlu * lung_plasma_dxorgluc - tr_lu_dxorgluc) / vlu_plasma
    d/dt(lung_dxorgluc) <- tr_lu_dxorgluc / vlu_tissue
    d/dt(rest_plasma_dxorgluc) <- (qre * arterial_dxorgluc - qre * rest_plasma_dxorgluc - tr_re_dxorgluc) / vre_plasma
    d/dt(rest_dxorgluc) <- tr_re_dxorgluc / vre_tissue
    d/dt(rest_vein_dxorgluc) <- (qre * rest_plasma_dxorgluc - qre * rest_vein_dxorgluc) / vrev
    d/dt(forearm_plasma_dxorgluc) <- (qfo * (1 - f_shunting_forearm) * arterial_dxorgluc - flow_fo_fov_dxorgluc - tr_fo_dxorgluc) / vfo_plasma
    d/dt(forearm_dxorgluc) <- tr_fo_dxorgluc / vfo_tissue
    d/dt(forearm_vein_dxorgluc) <- (flow_fo_fov_dxorgluc - qfo * forearm_vein_dxorgluc) / vfov
    d/dt(venous_dxorgluc) <- (qki * kidney_vein_dxorgluc + qh * hepatic_vein_dxorgluc + qre * rest_vein_dxorgluc + qfo * forearm_vein_dxorgluc - qlu * venous_dxorgluc) / vve
    d/dt(arterial_dxorgluc) <- (qlu * lung_plasma_dxorgluc - (qgu + qki + qha + qpa + qsp + qre + qfo * (1 - f_shunting_forearm)) * arterial_dxorgluc - flow_fo_fov_dxorgluc) / var_


    # ============ Observations ============
    # The paper samples plasma at the forearm (median cubital) vein. Reported
    # in nmol/L (= 1e6 * mmol/L concentration) for DXM, DXO and DXO-Glu.
    Cc <- 1e6 * forearm_vein
    Cc_dxor <- 1e6 * forearm_vein_dxor
    Cc_dxorgluc <- 1e6 * forearm_vein_dxorgluc

    # Cumulative urinary amounts (mmol) and the urinary cumulative metabolic
    # ratio UCMR = DXM / (DXO + DXO-Glu).
    Aurine_dxm <- urine
    Aurine_dxo <- urine_dxor
    Aurine_dxoglu <- urine_dxorgluc
    UCMR <- urine / (urine_dxor + urine_dxorgluc)

    Cc ~ prop(propSd)
  })
}
