Tensil_2026_losartan_pbpk <- function() {
  description <- paste(
    "PBPK/PD (whole-body, 49 ODEs, SBML/libroadrunner).",
    "Losartan, its active CYP2C9 metabolite E3174 (losartan carboxylic acid)",
    "and the E3174 glucuronide L158 in healthy adults, hepatic and renal",
    "impairment and CYP2C9 / ABCB1 genetic variants, built by Tensil,",
    "Myshkina and Koenig (2026) as a 'digital twin' from 25 curated clinical",
    "studies. The whole-body circulation links gut, liver, kidney and lung",
    "submodels through venous, arterial, portal and hepatic-venous plasma",
    "pools, with a lumped rest-of-body compartment. Oral losartan dissolves",
    "into the gut lumen; 60% of the absorbed flux enters the enterocytes",
    "(the rest is excreted in feces), from which P-glycoprotein (ABCB1)",
    "effluxes it back to the lumen or it passes to gut plasma. Liver tissue",
    "takes losartan up, converts it to E3174 by CYP2C9 (Michaelis-Menten)",
    "and E3174 to L158 (first order), exports E3174 and L158 to plasma and",
    "secretes all three species into bile, which empties into the gut lumen",
    "(enterohepatic circulation; E3174 and L158 then leave in feces). All",
    "three species are excreted in urine (first order from kidney plasma).",
    "Losartan alone partitions into lung and rest-of-body tissue. A",
    "renin-angiotensin-aldosterone (RAAS) layer driven by venous E3174",
    "raises renin secretion and suppresses aldosterone secretion through",
    "one Emax term; renin converts angiotensinogen to angiotensin I, which",
    "turns over to angiotensin II, and systolic, diastolic and mean",
    "arterial pressure move in proportion to the aldosterone deviation from",
    "baseline. Body weight, relative renal and liver function and CYP2C9",
    "activity are covariates; ABCB1 activity is a scaling parameter. The",
    "model is deterministic - the authors fit a typical individual and",
    "report no between-subject variability or residual error."
  )
  reference <- paste(
    "Tensil E, Myshkina M, Koenig M. A Digital Twin of the Angiotensin II",
    "Receptor Blocker Losartan: Physiologically Based Modeling of Blood",
    "Pressure Regulation. Pharmaceutics. 2026;18(2):262.",
    "doi:10.3390/pharmaceutics18020262.",
    "Optimized pharmacokinetic parameters from Supplementary Table S2,",
    "pharmacodynamic parameters from Supplementary Table S3 and the full",
    "parameter list (physiology, fixed literature constants, scaling",
    "factors) from Supplementary Table S7; fractional organ volumes, blood",
    "flows and the losartan partition coefficient from Methods Section 2.3.",
    "Rate laws and ODEs are not printed in the paper or supplement and were",
    "taken from the authors' archived SBML model, cited in Methods Section",
    "2.3 as the model source: Tensil E, Myshkina M, Koenig M.",
    "Physiologically Based Pharmacokinetic/Pharmacodynamic (PBPK/PD) Model",
    "of Losartan, v0.8.0. Zenodo. doi:10.5281/zenodo.18213663 (file",
    "models/losartan_body_flat.md, github.com/matthiaskoenig/losartan-model",
    "tag 0.8.0).",
    sep = " "
  )
  vignette <- "Tensil_2026_losartan_pbpk"
  units <- list(time = "min", dosing = "mg", concentration = "nmol/L")

  # Oral losartan plus intravenous losartan and intravenous E3174 (FDA 1995
  # study 60, Lo 1995), as simulated in the paper. The two intravenous dose
  # compartments are outside the auto-detected depot/central set, so they
  # are declared explicitly. An intravenous infusion is given as a `rate`
  # into the same compartments.
  dosing <- c("depot", "depot_iv", "depot_iv_e3174")

  # Drug states are CONCENTRATIONS in umol/L (amounts in umol for the
  # cumulative urine and feces pools); the RAAS hormones are in pmol/L. The
  # source SBML works in mmol and mmol/L throughout; the rescaling is exact
  # and keeps the solver away from its absolute tolerance (baseline renin is
  # 1e-9 mmol/L). Doses go into `depot` (oral, mg) or `depot_iv*`
  # (intravenous, mg) - never directly into a plasma or tissue compartment.
  paper_specific_compartments <- c(
    "depot_iv",
    "depot_iv_e3174",
    "gut_lumen_e3174",
    "gut_lumen_l158",
    "feces",
    "feces_e3174",
    "feces_l158",
    "gut_plasma",
    "gut_plasma_e3174",
    "gut_plasma_l158",
    "portal",
    "portal_e3174",
    "portal_l158",
    "liver_plasma",
    "liver_plasma_e3174",
    "liver_plasma_l158",
    "liver_e3174",
    "liver_l158",
    "hepatic_vein",
    "hepatic_vein_e3174",
    "hepatic_vein_l158",
    "kidney_plasma",
    "kidney_plasma_e3174",
    "kidney_plasma_l158",
    "urine_e3174",
    "urine_l158",
    "lung_plasma",
    "lung_plasma_e3174",
    "lung_plasma_l158",
    "rest",
    "rest_plasma",
    "rest_plasma_e3174",
    "rest_plasma_l158",
    "venous",
    "venous_e3174",
    "venous_l158",
    "arterial",
    "arterial_e3174",
    "arterial_l158",
    "renin",
    "ang1",
    "ang2",
    "aldosterone"
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Scales every absolute organ volume and, through cardiac output",
        "(COBW = 1.548 mL/s/kg), every blood flow: Methods Section 2.3",
        "states that 'absolute organ volumes and blood flows were",
        "calculated by scaling the corresponding fractional values with",
        "body weight'. Because the CYP2C9 Vmax and the hepatic and renal",
        "rate constants are expressed per litre of tissue, metabolic and",
        "excretory capacity scale with body weight as well. The reference",
        "individual is 75 kg (Supplementary Table S7, BW); the per-study",
        "simulations use the study mean bodyweight where reported",
        "(Supplementary Table S5)."
      ),
      source_name = "BW"
    ),
    RENALFUNC_REL = list(
      description = "Relative renal function as a fraction of normal",
      units = "(dimensionless)",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "1 = normal renal function. Multiplies the urinary excretion of",
        "losartan, E3174 and L158 (Methods Section 2.3: 'renal impairment",
        "was modeled as reduced clearance based on glomerular filtration",
        "rate following KDIGO guidelines'). The paper does not tabulate a",
        "KDIGO-to-scalar map; the Sica 1995 renal-impairment simulation",
        "uses 0.14, 0.5 and 0.95 for its three groups (Supplementary Table",
        "S5). Source name f_renal_function (SBML KI__f_renal_function)."
      ),
      source_name = "f_renal_function"
    ),
    HEPFUNC_REL = list(
      description = "Relative liver function as a fraction of normal",
      units = "(dimensionless)",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "1 = normal liver function. The source parameterises hepatic",
        "impairment with the OPPOSITE orientation, as a cirrhosis severity",
        "f_cirrhosis on [0, 0.95], so this column is the complement:",
        "f_cirrhosis = 1 - HEPFUNC_REL. f_cirrhosis acts twice, as the",
        "fraction of parenchymal liver volume lost (shrinking the",
        "metabolising tissue volume) and as the fraction of portal and",
        "hepatic-arterial blood shunted past the liver (Methods Section",
        "2.3, 'progressive cirrhosis'). The FDA 1995 study 67 cirrhosis",
        "simulation uses f_cirrhosis = 0.67 (HEPFUNC_REL = 0.33;",
        "Supplementary Table S5) and the paper's Figure 4 scans",
        "f_cirrhosis over 0-0.9. Take care with the inversion when",
        "reproducing the paper's figures, which are labelled by",
        "f_cirrhosis."
      ),
      source_name = "f_cirrhosis"
    ),
    CYP2C9 = list(
      description = "CYP2C9 activity relative to the *1/*1 genotype",
      units = "(fraction of *1/*1 activity)",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "1 = *1/*1 (normal) activity. Multiplies the hepatic CYP2C9",
        "losartan -> E3174 Vmax directly; no centering is applied (Methods",
        "Section 2.3, 'allele-specific activity scaling'). The per-study",
        "simulations use 1.0 for *1/*1 and, for reduced-function carriers,",
        "0.8, 0.6, 0.585, 0.525, 0.385 and 0.17 (Supplementary Table S5:",
        "Bae 2011, Han 2009, Huang 2021, Lee 2003, Li 2009, Sekino 2003,",
        "Yasar 2002). The paper does not print the genotype-to-activity",
        "map, so which value belongs to which diplotype is taken from the",
        "study groups in the vignette. Source name f_cyp2c9 (SBML",
        "LI__f_cyp2c9)."
      ),
      source_name = "f_cyp2c9"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "losartan", units = "mg", specimen = "administration site", verified = TRUE),
    depot_iv = list(analyte = "losartan", units = "mg", specimen = "administration site", verified = TRUE),
    depot_iv_e3174 = list(analyte = "E3174", units = "mg", specimen = "administration site", verified = TRUE),
    gut_lumen = list(analyte = "losartan", units = "umol/L", specimen = "administration site", verified = TRUE),
    gut_lumen_e3174 = list(analyte = "E3174", units = "umol/L", specimen = "administration site", verified = TRUE),
    gut_lumen_l158 = list(analyte = "L158", units = "umol/L", specimen = "administration site", verified = TRUE),
    gut = list(analyte = "losartan", units = "umol/L", specimen = "tissue", verified = TRUE),
    feces = list(analyte = "losartan", units = "umol", specimen = "faeces", verified = TRUE),
    feces_e3174 = list(analyte = "E3174", units = "umol", specimen = "faeces", verified = TRUE),
    feces_l158 = list(analyte = "L158", units = "umol", specimen = "faeces", verified = TRUE),
    gut_plasma = list(analyte = "losartan", units = "umol/L", specimen = "plasma", verified = TRUE),
    gut_plasma_e3174 = list(analyte = "E3174", units = "umol/L", specimen = "plasma", verified = TRUE),
    gut_plasma_l158 = list(analyte = "L158", units = "umol/L", specimen = "plasma", verified = TRUE),
    portal = list(analyte = "losartan", units = "umol/L", specimen = "plasma", verified = TRUE),
    portal_e3174 = list(analyte = "E3174", units = "umol/L", specimen = "plasma", verified = TRUE),
    portal_l158 = list(analyte = "L158", units = "umol/L", specimen = "plasma", verified = TRUE),
    liver_plasma = list(analyte = "losartan", units = "umol/L", specimen = "plasma", verified = TRUE),
    liver_plasma_e3174 = list(analyte = "E3174", units = "umol/L", specimen = "plasma", verified = TRUE),
    liver_plasma_l158 = list(analyte = "L158", units = "umol/L", specimen = "plasma", verified = TRUE),
    liver = list(analyte = "losartan", units = "umol/L", specimen = "tissue", verified = TRUE),
    liver_e3174 = list(analyte = "E3174", units = "umol/L", specimen = "tissue", verified = TRUE),
    liver_l158 = list(analyte = "L158", units = "umol/L", specimen = "tissue", verified = TRUE),
    hepatic_vein = list(analyte = "losartan", units = "umol/L", specimen = "plasma", verified = TRUE),
    hepatic_vein_e3174 = list(analyte = "E3174", units = "umol/L", specimen = "plasma", verified = TRUE),
    hepatic_vein_l158 = list(analyte = "L158", units = "umol/L", specimen = "plasma", verified = TRUE),
    kidney_plasma = list(analyte = "losartan", units = "umol/L", specimen = "plasma", verified = TRUE),
    kidney_plasma_e3174 = list(analyte = "E3174", units = "umol/L", specimen = "plasma", verified = TRUE),
    kidney_plasma_l158 = list(analyte = "L158", units = "umol/L", specimen = "plasma", verified = TRUE),
    urine = list(analyte = "losartan", units = "umol", specimen = "urine", verified = TRUE),
    urine_e3174 = list(analyte = "E3174", units = "umol", specimen = "urine", verified = TRUE),
    urine_l158 = list(analyte = "L158", units = "umol", specimen = "urine", verified = TRUE),
    lung_plasma = list(analyte = "losartan", units = "umol/L", specimen = "plasma", verified = TRUE),
    lung_plasma_e3174 = list(analyte = "E3174", units = "umol/L", specimen = "plasma", verified = TRUE),
    lung_plasma_l158 = list(analyte = "L158", units = "umol/L", specimen = "plasma", verified = TRUE),
    lung = list(analyte = "losartan", units = "umol/L", specimen = "tissue", verified = TRUE),
    rest_plasma = list(analyte = "losartan", units = "umol/L", specimen = "plasma", verified = TRUE),
    rest_plasma_e3174 = list(analyte = "E3174", units = "umol/L", specimen = "plasma", verified = TRUE),
    rest_plasma_l158 = list(analyte = "L158", units = "umol/L", specimen = "plasma", verified = TRUE),
    rest = list(analyte = "losartan", units = "umol/L", specimen = "tissue", verified = TRUE),
    venous = list(analyte = "losartan", units = "umol/L", specimen = "plasma", verified = TRUE),
    venous_e3174 = list(analyte = "E3174", units = "umol/L", specimen = "plasma", verified = TRUE),
    venous_l158 = list(analyte = "L158", units = "umol/L", specimen = "plasma", verified = TRUE),
    arterial = list(analyte = "losartan", units = "umol/L", specimen = "plasma", verified = TRUE),
    arterial_e3174 = list(analyte = "E3174", units = "umol/L", specimen = "plasma", verified = TRUE),
    arterial_l158 = list(analyte = "L158", units = "umol/L", specimen = "plasma", verified = TRUE),
    renin = list(analyte = "renin", units = "pmol/L", specimen = "plasma", verified = TRUE),
    ang1 = list(analyte = "angiotensin I", units = "pmol/L", specimen = "plasma", verified = TRUE),
    ang2 = list(analyte = "angiotensin II", units = "pmol/L", specimen = "plasma", verified = TRUE),
    aldosterone = list(analyte = "aldosterone", units = "pmol/L", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = NA_integer_,
    n_studies = 25L,
    age_range = "adults; not tabulated per study in the source",
    weight_range = paste(
      "study-specific mean bodyweight where reported (52-84.6 kg across",
      "Supplementary Table S5); 75 kg reference individual otherwise"
    ),
    sex_female_pct = NA_real_,
    disease_state = paste(
      "Healthy volunteers, CYP2C9 and ABCB1 genotyped volunteers, and",
      "subjects with hepatic (cirrhosis) or renal impairment (Table 1).",
      "Animal, pediatric, single-patient, cocktail and combination studies",
      "were excluded (Methods Section 2.1)."
    ),
    dose_range = paste(
      "Oral losartan 2-200 mg single and multiple dose (Table 1);",
      "intravenous losartan 10-30 mg and intravenous E3174 10-20 mg",
      "(FDA 1995 studies 60 and 67, Lo 1995)."
    ),
    regions = "not reported; the curated studies are international",
    notes = paste(
      "Deterministic typical-individual model. Methods Section 2.3 states",
      "that 'all simulations were performed deterministically using the",
      "optimized parameter set representing the typical (mean) individual'",
      "and that inter-individual variability was NOT included, so no IIV",
      "and no residual-error model are reported and none are invented",
      "here. Data were curated into PK-DB (identifiers in Table 1) from a",
      "PubMed and PKPDAI search screening 787 records down to 25 studies.",
      "Parameters were fitted sequentially by weighted least squares with",
      "100 multi-start local optimisation runs - pharmacokinetic first",
      "(12 parameters, Table S2), then pharmacodynamic (4 parameters,",
      "Table S3). Pharmacodynamic data (renin, angiotensin, aldosterone,",
      "blood pressure) came from a smaller subset of studies and showed",
      "high between-study variability (Discussion)."
    )
  )

  ini({
    # ================= Reference physiology ==============================
    # Supplementary Table S7 and Methods Section 2.3 give the fractional
    # volumes and blood flows; the blood-pool fractions FVhv / FVpo, Fblood,
    # the enterocyte volume and the RAAS distribution volume are not printed
    # and come from the source SBML.
    covbw    <- fixed(1.548); label("Cardiac output per unit body weight (mL/s/kg)")  # Table S7 COBW = 1.55; SBML COBW = 1.548
    fcardio  <- fixed(1);     label("Cardiac function scaling factor (unitless)")     # SBML f_cardiac_function, normal = 1
    hct      <- fixed(0.51);  label("Hematocrit (fraction)")                          # Table S7 HCT = 0.51
    fblood   <- fixed(0.02);  label("Blood fraction of organ volume (unitless)")      # SBML Fblood
    v_entero <- fixed(1);     label("Enterocyte (intestinal lining) volume (L)")      # SBML GU__Ventero
    vplasma  <- fixed(5);     label("Plasma volume of the RAAS hormones (L)")         # SBML Vplasma

    fvar <- fixed(0.0257); label("Arterial blood fractional volume (L/kg)")  # Table S7 FVar = 0.0257
    fvgu <- fixed(0.0171); label("Gut fractional volume (L/kg)")             # Methods 2.3 'FVgu=1.71%'; Table S7 FVgu
    fvhv <- fixed(0.001);  label("Hepatic venous fractional volume (L/kg)")  # SBML FVhv
    fvki <- fixed(0.0044); label("Kidney fractional volume (L/kg)")          # Methods 2.3 'FVki=0.44%'; Table S7 FVki
    fvli <- fixed(0.021);  label("Liver fractional volume (L/kg)")           # Methods 2.3 'FVli=2.10%'; Table S7 FVli
    fvlu <- fixed(0.0076); label("Lung fractional volume (L/kg)")            # Methods 2.3 'FVlu=0.76%'; Table S7 FVlu
    fvpo <- fixed(0.001);  label("Portal fractional volume (L/kg)")          # SBML FVpo
    fvve <- fixed(0.0514); label("Venous blood fractional volume (L/kg)")    # Table S7 FVve = 0.0514

    fqgu <- fixed(0.18);  label("Gut fractional blood flow (unitless)")                  # Methods 2.3 'FQgu=18.00%'; Table S7 FQgu
    fqh  <- fixed(0.215); label("Hepatic venous-side fractional blood flow (unitless)")  # Methods 2.3 'FQh=21.50%'; Table S7 FQh
    fqki <- fixed(0.19);  label("Kidney fractional blood flow (unitless)")               # Methods 2.3 'FQki=19.00%'; Table S7 FQki
    fqlu <- fixed(1);     label("Lung fractional blood flow (unitless)")                 # Methods 2.3 'FQlu=100%'

    # ================= Molecular weights =================================
    mr_los   <- fixed(422.911); label("Molecular weight of losartan (g/mol)")  # SBML Mr_los
    mr_e3174 <- fixed(436.894); label("Molecular weight of E3174 (g/mol)")     # SBML Mr_e3174

    # ================= Distribution (losartan only) ======================
    # Methods Section 2.3: one partition coefficient for all tissues for
    # losartan; no tissue partitioning for E3174 and L158.
    kp_los      <- fixed(3.26214895692923);  label("Tissue-to-plasma partition coefficient, losartan (unitless)")  # Table S2 Kp_los = 3.262; Methods 2.3 'KpLOS=3.262'
    ftissue_los <- fixed(0.144942664253079); label("Tissue distribution rate, losartan (L/min)")                   # Table S2 ftissue_los = 0.145

    # ================= Absorption ========================================
    ka_dis_los <- fixed(2);                  label("Dissolution rate of losartan in the GI tract (1/hr)")               # Table S7 GU__Ka_dis_los = 2
    losabs_k   <- fixed(0.0361772268991116); label("Absorption rate of losartan in the GI tract (1/min)")                # Table S2 GU__LOSABS_k = 0.036
    f_losefl_k <- fixed(1.02003554075262);   label("P-glycoprotein efflux rate relative to absorption (unitless)")      # Table S2 GU__f_LOSEFL_k = 1.020
    f_los_abs  <- fixed(0.6);                label("Fraction of the luminal absorption flux entering enterocytes (unitless)")  # Table S7 GU__F_los_abs = 0.6
    metexc_k   <- fixed(2.69641608199351e-05); label("Fecal excretion rate of E3174 and L158 from the gut lumen (1/min)")  # Table S2 GU__METEXC_k = 2.696e-5
    f_oatp2b1  <- fixed(1);                  label("Intestinal uptake (OATP2B1) activity scaling factor (unitless)")    # Table S7 GU__f_OATP2B1 = 1
    f_abcb1    <- fixed(1);                  label("Intestinal P-glycoprotein (ABCB1) activity scaling factor (unitless)")  # Table S7 GU__f_abcb1 = 1

    # ================= Liver =============================================
    losim_k         <- fixed(1);                   label("Hepatic losartan uptake rate (1/min)")                         # SBML LI__LOSIM_k (not fitted, not tabulated)
    los2e3174_vmax  <- fixed(0.00072514631731673); label("Vmax, CYP2C9 conversion of losartan to E3174 (mmol/min/L)")   # Table S2 LI__LOS2E3174_Vmax = 7.251e-4
    los2e3174_km    <- fixed(0.0011);              label("Km, CYP2C9 conversion of losartan to E3174 (mmol/L)")          # Table S7 LI__LOS2E3174_Km_los = 0.0011
    e3174l158_k     <- fixed(0.00113400909936956); label("Hepatic conversion rate of E3174 to L158 (1/min)")             # Table S2 LI__E3174L158_k = 1.134e-3
    li_e3174ex_k    <- fixed(0.0112001182272967);  label("Hepatic export rate of E3174 to plasma (1/min)")               # Table S2 LI__E3174EX_k = 0.011
    li_l158ex_k     <- fixed(1);                   label("Hepatic export rate of L158 to plasma (1/min)")                # SBML LI__L158EX_k (not fitted, not tabulated)
    mbiex_k         <- fixed(0.066315021777011);   label("Biliary export rate of losartan, E3174 and L158 (1/min)")       # Table S2 LI__MBIEX_k = 0.066

    # ================= Kidney ============================================
    ki_losex_k   <- fixed(0.0773896174258953); label("Urinary excretion rate of losartan (1/min)")  # Table S2 KI__LOSEX_k = 0.077
    ki_e3174ex_k <- fixed(0.0289271732633318); label("Urinary excretion rate of E3174 (1/min)")     # Table S2 KI__E3174EX_k = 0.029
    ki_l158ex_k  <- fixed(0.28890834506773);   label("Urinary excretion rate of L158 (1/min)")      # Table S2 KI__L158EX_k = 0.289

    # ================= Pharmacodynamics (RAAS) ===========================
    anggen2ang1_k   <- fixed(0.100300143544081);  label("Angiotensinogen to angiotensin I conversion rate (L/min)")  # Table S3 ANGGEN2ANG1_k = 0.100
    e50_e3174       <- fixed(0.000291072740111984); label("Venous E3174 giving half-maximal RAAS effect (mmol/L)")   # Table S3 E50_e3174 = 2.911e-4
    aldsec_k        <- fixed(1.0105027952686e-06); label("Aldosterone secretion rate at baseline (mmol/min)")       # Table S7 ALDSEC_k = 1.01e-06 (Table S3 prints 1.011 - see vignette Errata)
    bp_ald_fe       <- fixed(0.311999262694736);  label("Effect of aldosterone on blood pressure (unitless)")       # Table S3 BP_ald_fe = 0.312
    rensec_k        <- fixed(0.1);                label("Renin secretion rate at baseline (mmol/min)")              # Table S7 RENSEC_k = 0.1
    rensec_fa_e3174 <- fixed(5);                  label("Maximal fold activation of renin secretion by E3174 (unitless)")  # Table S7 RENSEC_fa_e3174 = 5

    # Baseline (reference) hormone levels and blood pressures. The per-study
    # simulations reset these to each study's placebo/baseline values
    # (Supplementary Table S5); the source SBML stores them in mmol/L, here
    # they are in pmol/L (SBML value x 1e9).
    ren_ref    <- fixed(1);   label("Baseline plasma renin (pmol/L)")           # SBML ren_ref = 1e-9 mmol/L
    ang1_ref   <- fixed(10);  label("Baseline plasma angiotensin I (pmol/L)")   # SBML ang1_ref = 1e-8 mmol/L
    ang2_ref   <- fixed(8);   label("Baseline plasma angiotensin II (pmol/L)")  # SBML ang2_ref = 8e-9 mmol/L
    ald_ref    <- fixed(250); label("Baseline plasma aldosterone (pmol/L)")     # SBML ald_ref = 2.5e-7 mmol/L
    anggen_ref <- fixed(10);  label("Plasma angiotensinogen, held constant (pmol/L)")  # SBML anggen_ref = 1e-8 mmol/L
    sbp_ref    <- fixed(120); label("Baseline systolic blood pressure (mmHg)")  # Results 3.2 'physiological reference values of 120/80 mmHg'
    dbp_ref    <- fixed(80);  label("Baseline diastolic blood pressure (mmHg)") # Results 3.2 'physiological reference values of 120/80 mmHg'

    # ================= Intravenous input =================================
    ti_los   <- fixed(10); label("Nominal intravenous injection time, losartan (s)")  # SBML ti_los; first-order release with half-life ti_los
    ti_e3174 <- fixed(10); label("Nominal intravenous injection time, E3174 (s)")     # SBML ti_e3174

    # ================= Residual error ====================================
    # Methods Section 2.3: deterministic typical-individual simulations, no
    # population variability; no residual-error model is reported, so it is
    # fixed at zero rather than invented.
    propSd <- fixed(0); label("Proportional residual error (fraction; ZERO - not reported by the source)")  # Methods 2.3: deterministic simulations
  })

  model({
    # ================= Unit conversions ==================================
    # Drug states in umol/L (source: mmol/L) and RAAS states in pmol/L
    # (source: mmol/L).
    cf_um <- 1000
    cf_pm <- 1e9

    # ================= Scenario scalars ==================================
    f_renal       <- RENALFUNC_REL
    f_cyp2c9      <- CYP2C9
    f_cirrhosis   <- 1 - HEPFUNC_REL
    f_shunts      <- f_cirrhosis
    f_tissue_loss <- f_cirrhosis

    # ================= Volumes ===========================================
    fvre <- 1 - (fvgu + fvki + fvli + fvlu + fvve + fvar)

    vgu <- WT * fvgu
    vki <- WT * fvki
    vli <- WT * fvli
    vlu <- WT * fvlu
    vre <- WT * fvre

    # Arterial and venous blood pools are corrected for the blood held within
    # the organ volumes.
    var_ <- WT * fvar - (fvar / (fvar + fvve)) * WT * fblood * (1 - fvve - fvar)
    vve  <- WT * fvve - (fvve / (fvar + fvve)) * WT * fblood * (1 - fvve - fvar)
    vpo  <- (1 - hct) * (WT * fvpo - (fvpo / (fvar + fvve + fvpo + fvhv)) * WT * fblood * (1 - (fvar + fvve + fvpo + fvhv)))
    vhv  <- (1 - hct) * (WT * fvhv - (fvhv / (fvar + fvve + fvpo + fvhv)) * WT * fblood * (1 - (fvar + fvve + fvpo + fvhv)))

    vgu_plasma <- vgu * fblood * (1 - hct)
    vki_plasma <- vki * fblood * (1 - hct)
    vli_plasma <- vli * fblood * (1 - hct)
    vlu_plasma <- vlu * fblood * (1 - hct)
    vre_plasma <- vre * fblood * (1 - hct)

    vki_tissue <- vki * (1 - fblood)
    vlu_tissue <- vlu * (1 - fblood)
    vre_tissue <- vre * (1 - fblood)
    # Cirrhosis destroys functional parenchyma, so the metabolising liver
    # tissue volume shrinks with severity.
    vli_tissue <- vli * (1 - f_tissue_loss) * (1 - fblood)

    # ================= Blood flows =======================================
    co <- fcardio * WT * covbw
    qc <- (co / 1000) * 60

    fqre <- 1 - (fqki + fqh)

    qgu <- qc * fqgu
    qh  <- qc * fqh
    qki <- qc * fqki
    qlu <- qc * fqlu
    qre <- qc * fqre
    qha <- qh - qgu
    qpo <- qgu

    # ================= Intestine =========================================
    # Rates in umol/min. Losartan leaving the lumen is split: f_los_abs
    # (scaled by the uptake activity) enters the enterocytes, the rest is
    # excreted in feces. P-gp efflux and export to gut plasma are both
    # written on the enterocyte concentration times the GUT volume, as in
    # the source.
    dissolution <- cf_um * (ka_dis_los / 60) * depot / mr_los
    absorption  <- losabs_k * vgu * gut_lumen
    losabs      <- f_oatp2b1 * f_los_abs * absorption
    losexc      <- (1 - f_los_abs) * absorption
    losefl      <- f_abcb1 * f_losefl_k * losabs_k * vgu * gut
    losex       <- losabs_k * vgu * gut
    e3174exc    <- metexc_k * vgu * gut_lumen_e3174
    l158exc     <- metexc_k * vgu * gut_lumen_l158

    # ================= Intravenous input =================================
    ki_los   <- (0.693 / ti_los) * 60
    ki_e3174 <- (0.693 / ti_e3174) * 60
    iv_los   <- cf_um * ki_los * depot_iv / mr_los
    iv_e3174 <- cf_um * ki_e3174 * depot_iv_e3174 / mr_e3174

    # ================= Liver =============================================
    losim      <- losim_k * vli_tissue * (liver_plasma - liver)
    los2e3174  <- f_cyp2c9 * cf_um * los2e3174_vmax * vli_tissue * liver /
      (liver + cf_um * los2e3174_km)
    e3174l158  <- e3174l158_k * vli_tissue * liver_e3174
    li_e3174ex <- li_e3174ex_k * vli_tissue * (liver_e3174 - liver_plasma_e3174)
    li_l158ex  <- li_l158ex_k * vli_tissue * (liver_l158 - liver_plasma_l158)
    # Bile is a pass-through: everything secreted into bile reaches the gut
    # lumen at the same instant (the source's bile pools have a zero net
    # rate and are omitted).
    losbiex    <- mbiex_k * vli_tissue * liver
    e3174biex  <- mbiex_k * vli_tissue * liver_e3174
    l158biex   <- mbiex_k * vli_tissue * liver_l158

    # ================= Kidney ============================================
    ki_losex   <- f_renal * ki_losex_k * vki_tissue * kidney_plasma
    ki_e3174ex <- f_renal * ki_e3174ex_k * vki_tissue * kidney_plasma_e3174
    ki_l158ex  <- f_renal * ki_l158ex_k * vki_tissue * kidney_plasma_l158

    # ================= Perfusion-limited tissue distribution =============
    transport_lu <- ftissue_los * (lung_plasma * kp_los - lung)
    transport_re <- ftissue_los * (rest_plasma * kp_los - rest)

    # ================= Pharmacodynamics (RAAS) ===========================
    # One Emax term on venous E3174 (both in mmol/L) raises renin secretion
    # and suppresses aldosterone secretion. Every degradation constant is
    # set so that the hormones sit at their reference levels without drug.
    fe_e3174 <- (venous_e3174 / cf_um) / (e50_e3174 + venous_e3174 / cf_um)

    rensec      <- cf_pm * rensec_k * (1 + rensec_fa_e3174 * fe_e3174)
    rendeg      <- cf_pm * rensec_k / ren_ref * renin
    anggen2ang1 <- anggen2ang1_k * anggen_ref * renin / ren_ref
    ang1ang2    <- anggen2ang1_k * anggen_ref / ang1_ref * ang1
    ang2deg     <- anggen2ang1_k * anggen_ref / ang2_ref * ang2
    aldsec      <- cf_pm * aldsec_k * (1 - fe_e3174)
    alddeg      <- cf_pm * aldsec_k / ald_ref * aldosterone

    # ================= ODEs ==============================================
    # Dosing compartments (mg).
    d/dt(depot)          <- -dissolution * mr_los / cf_um
    d/dt(depot_iv)       <- -iv_los * mr_los / cf_um
    d/dt(depot_iv_e3174) <- -iv_e3174 * mr_e3174 / cf_um

    # Gut lumen receives dissolved drug and bile; enterocytes.
    d/dt(gut_lumen)       <- (dissolution + losbiex - losabs + losefl - losexc) / vgu
    d/dt(gut_lumen_e3174) <- (e3174biex - e3174exc) / vgu
    d/dt(gut_lumen_l158)  <- (l158biex - l158exc) / vgu
    d/dt(gut)             <- (losabs - losefl - losex) / v_entero

    d/dt(feces)       <- losexc
    d/dt(feces_e3174) <- e3174exc
    d/dt(feces_l158)  <- l158exc

    # Gut plasma drains to the portal vein.
    d/dt(gut_plasma)       <- (qgu * arterial - qgu * gut_plasma + losex) / vgu_plasma
    d/dt(gut_plasma_e3174) <- (qgu * arterial_e3174 - qgu * gut_plasma_e3174) / vgu_plasma
    d/dt(gut_plasma_l158)  <- (qgu * arterial_l158 - qgu * gut_plasma_l158) / vgu_plasma

    # Portal vein; cirrhotic shunts divert part of it into the hepatic vein.
    d/dt(portal)       <- (qgu * gut_plasma - qpo * portal) / vpo
    d/dt(portal_e3174) <- (qgu * gut_plasma_e3174 - qpo * portal_e3174) / vpo
    d/dt(portal_l158)  <- (qgu * gut_plasma_l158 - qpo * portal_l158) / vpo

    # Liver plasma.
    d/dt(liver_plasma) <- ((1 - f_shunts) * qha * arterial + (1 - f_shunts) * qpo * portal -
      (1 - f_shunts) * (qpo + qha) * liver_plasma - losim) / vli_plasma
    d/dt(liver_plasma_e3174) <- ((1 - f_shunts) * qha * arterial_e3174 +
      (1 - f_shunts) * qpo * portal_e3174 -
      (1 - f_shunts) * (qpo + qha) * liver_plasma_e3174 + li_e3174ex) / vli_plasma
    d/dt(liver_plasma_l158) <- ((1 - f_shunts) * qha * arterial_l158 +
      (1 - f_shunts) * qpo * portal_l158 -
      (1 - f_shunts) * (qpo + qha) * liver_plasma_l158 + li_l158ex) / vli_plasma

    # Liver tissue.
    d/dt(liver)       <- (losim - los2e3174 - losbiex) / vli_tissue
    d/dt(liver_e3174) <- (los2e3174 - e3174l158 - li_e3174ex - e3174biex) / vli_tissue
    d/dt(liver_l158)  <- (e3174l158 - li_l158ex - l158biex) / vli_tissue

    # Hepatic vein collects liver outflow and both shunts.
    d/dt(hepatic_vein) <- (f_shunts * qha * arterial + f_shunts * qpo * portal +
      (1 - f_shunts) * (qpo + qha) * liver_plasma - qh * hepatic_vein) / vhv
    d/dt(hepatic_vein_e3174) <- (f_shunts * qha * arterial_e3174 +
      f_shunts * qpo * portal_e3174 +
      (1 - f_shunts) * (qpo + qha) * liver_plasma_e3174 - qh * hepatic_vein_e3174) / vhv
    d/dt(hepatic_vein_l158) <- (f_shunts * qha * arterial_l158 +
      f_shunts * qpo * portal_l158 +
      (1 - f_shunts) * (qpo + qha) * liver_plasma_l158 - qh * hepatic_vein_l158) / vhv

    # Kidney plasma and urine.
    d/dt(kidney_plasma)       <- (qki * arterial - qki * kidney_plasma - ki_losex) / vki_plasma
    d/dt(kidney_plasma_e3174) <- (qki * arterial_e3174 - qki * kidney_plasma_e3174 - ki_e3174ex) / vki_plasma
    d/dt(kidney_plasma_l158)  <- (qki * arterial_l158 - qki * kidney_plasma_l158 - ki_l158ex) / vki_plasma

    d/dt(urine)       <- ki_losex
    d/dt(urine_e3174) <- ki_e3174ex
    d/dt(urine_l158)  <- ki_l158ex

    # Lung sits between the venous and arterial pools.
    d/dt(lung_plasma)       <- (qlu * venous - qlu * lung_plasma - transport_lu) / vlu_plasma
    d/dt(lung_plasma_e3174) <- (qlu * venous_e3174 - qlu * lung_plasma_e3174) / vlu_plasma
    d/dt(lung_plasma_l158)  <- (qlu * venous_l158 - qlu * lung_plasma_l158) / vlu_plasma
    d/dt(lung)              <- transport_lu / vlu_tissue

    # Rest of body.
    d/dt(rest_plasma)       <- (qre * arterial - qre * rest_plasma - transport_re) / vre_plasma
    d/dt(rest_plasma_e3174) <- (qre * arterial_e3174 - qre * rest_plasma_e3174) / vre_plasma
    d/dt(rest_plasma_l158)  <- (qre * arterial_l158 - qre * rest_plasma_l158) / vre_plasma
    d/dt(rest)              <- transport_re / vre_tissue

    # Systemic circulation.
    d/dt(venous) <- (iv_los + qki * kidney_plasma + qh * hepatic_vein +
      qre * rest_plasma - qlu * venous) / vve
    d/dt(venous_e3174) <- (iv_e3174 + qki * kidney_plasma_e3174 + qh * hepatic_vein_e3174 +
      qre * rest_plasma_e3174 - qlu * venous_e3174) / vve
    d/dt(venous_l158) <- (qki * kidney_plasma_l158 + qh * hepatic_vein_l158 +
      qre * rest_plasma_l158 - qlu * venous_l158) / vve

    d/dt(arterial) <- (qlu * lung_plasma - qki * arterial - qha * arterial -
      qgu * arterial - qre * arterial) / var_
    d/dt(arterial_e3174) <- (qlu * lung_plasma_e3174 - qki * arterial_e3174 -
      qha * arterial_e3174 - qgu * arterial_e3174 - qre * arterial_e3174) / var_
    d/dt(arterial_l158) <- (qlu * lung_plasma_l158 - qki * arterial_l158 -
      qha * arterial_l158 - qgu * arterial_l158 - qre * arterial_l158) / var_

    # RAAS hormones (pmol/L), starting at their reference levels.
    d/dt(renin)       <- (rensec - rendeg) / vplasma
    d/dt(ang1)        <- (anggen2ang1 - ang1ang2) / vplasma
    d/dt(ang2)        <- (ang1ang2 - ang2deg) / vplasma
    d/dt(aldosterone) <- (aldsec - alddeg) / vplasma
    renin(0)       <- ren_ref
    ang1(0)        <- ang1_ref
    ang2(0)        <- ang2_ref
    aldosterone(0) <- ald_ref

    # ================= Observations ======================================
    # Venous plasma concentrations in nmol/L, matching the paper's figures.
    Cc       <- 1000 * venous
    Cc_e3174 <- 1000 * venous_e3174
    Cc_l158  <- 1000 * venous_l158

    # Blood pressure moves in proportion to the aldosterone deviation from
    # baseline (Results 3.2 and Discussion).
    SBP <- sbp_ref + bp_ald_fe * sbp_ref * (aldosterone - ald_ref) / ald_ref
    DBP <- dbp_ref + bp_ald_fe * dbp_ref * (aldosterone - ald_ref) / ald_ref
    MAP <- DBP + (SBP - DBP) / 3

    Cc ~ prop(propSd)
  })
}
