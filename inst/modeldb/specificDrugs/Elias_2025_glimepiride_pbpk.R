Elias_2025_glimepiride_pbpk <- function() {
  description <- paste(
    "PBPK (whole-body, 41 ODEs, SBML/libroadrunner).",
    "Glimepiride and its two CYP2C9-pathway metabolites M1",
    "(cyclohexylhydroxymethyl-glimepiride) and M2",
    "(carboxy-glimepiride) in healthy adults, type 2 diabetes mellitus and",
    "renal or hepatic impairment, built by Elias and Koenig (2025) as a",
    "'digital twin' from 20 curated clinical studies. The whole-body",
    "circulation links gut, liver, kidney and lung submodels through venous,",
    "arterial, portal and hepatic-venous plasma pools, with a lumped",
    "rest-of-body compartment. Oral drug dissolves from the dose compartment",
    "into the gut lumen and is absorbed first-order into gut plasma.",
    "Glimepiride is taken up into liver tissue and converted to M1 by",
    "CYP2C9 (Michaelis-Menten); M1 is converted to M2 (first order) and",
    "both metabolites are exported back to plasma, excreted in urine",
    "(first order from kidney plasma) and secreted into the gut lumen for",
    "fecal elimination. Unchanged glimepiride is not excreted. Tissue",
    "distribution into lung and rest-of-body uses one partition coefficient",
    "and one distribution rate shared by all three species. Four scenario",
    "covariates reproduce the paper's variability analyses: body weight",
    "(scales every organ volume and blood flow), relative renal function",
    "(scales metabolite urinary excretion), relative liver function",
    "(cirrhotic parenchymal loss and portal/arterial shunting) and CYP2C9",
    "activity relative to *1/*1 (scales the M1-formation Vmax). The model is",
    "deterministic - the authors fit a typical individual and report no",
    "between-subject variability or residual error."
  )
  reference <- paste(
    "Elias M, Koenig M. A digital twin of glimepiride for personalized and",
    "stratified diabetes treatment. Front Pharmacol. 2025;16:1686415.",
    "doi:10.3389/fphar.2025.1686415.",
    "Optimized parameters from Supplementary Table S1; submodel rate laws",
    "from Supplementary Section S3 (Intestine, Kidney and Liver Models);",
    "renal, hepatic and CYP2C9 scaling from Methods Section 4.2 and",
    "Supplementary Table S2. The whole-body circulation ODEs, physiology",
    "(fractional organ volumes and blood flows), the fixed Michaelis",
    "constant, molecular weights and dissolution rate are not printed in",
    "the paper or supplement and were taken from the authors' archived SBML",
    "model, cited in Methods Section 4.2 and Supplementary S3 as the",
    "definitive model source: Elias M, Koenig M. Physiologically based",
    "pharmacokinetic (PBPK) model of glimepiride, v0.6.1. Zenodo.",
    "doi:10.5281/zenodo.15189579 (file models/glimepiride_body_flat.md,",
    "github.com/matthiaskoenig/glimepiride-model tag 0.6.1).",
    sep = " "
  )
  vignette <- "Elias_2025_glimepiride_pbpk"
  units <- list(time = "min", dosing = "mg", concentration = "umol/L")

  # Oral glimepiride plus intravenous glimepiride (Badian 1994, FDA 1995,
  # Ratheiser 1993) and intravenous M1 (Badian 1996), as simulated in the
  # paper; intravenous M2 exists in the source SBML and is carried for
  # completeness.
  dosing <- c("depot", "depot_iv", "depot_iv_m1", "depot_iv_m2")

  # States are CONCENTRATIONS (mmol/L) throughout the circulation and the
  # tissues, following the source SBML; the dosing compartments are amounts
  # in mg and the urine/feces pools are cumulative amounts in mmol. Doses go
  # into `depot` (oral) or `depot_iv*` (intravenous) - never directly into a
  # plasma or tissue compartment.
  paper_specific_compartments <- c(
    "depot_iv",
    "depot_iv_m1",
    "depot_iv_m2",
    "gut_lumen",
    "gut_lumen_m1",
    "gut_lumen_m2",
    "feces_m1",
    "feces_m2",
    "gut_plasma",
    "gut_plasma_m1",
    "gut_plasma_m2",
    "portal",
    "portal_m1",
    "portal_m2",
    "liver_plasma",
    "liver_plasma_m1",
    "liver_plasma_m2",
    "liver_m1",
    "liver_m2",
    "hepatic_vein",
    "hepatic_vein_m1",
    "hepatic_vein_m2",
    "kidney_plasma",
    "kidney_plasma_m1",
    "kidney_plasma_m2",
    "urine_m1",
    "urine_m2",
    "lung_plasma",
    "lung_plasma_m1",
    "lung_plasma_m2",
    "lung_m1",
    "lung_m2",
    "rest",
    "rest_m1",
    "rest_m2",
    "rest_plasma",
    "rest_plasma_m1",
    "rest_plasma_m2",
    "venous",
    "venous_m1",
    "venous_m2",
    "arterial",
    "arterial_m1",
    "arterial_m2"
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Scales every absolute organ volume and, through cardiac output",
        "(COBW = 1.548 mL/s/kg), every blood flow; because the CYP2C9 Vmax",
        "is expressed per litre of liver tissue, the metabolic capacity",
        "scales with body weight as well (Methods Section 4.2: 'Bodyweight",
        "effects were incorporated by scaling organ volumes, blood flows,",
        "and metabolic rates'). The reference individual is 75 kg (source",
        "SBML parameter BW)."
      ),
      source_name = "BW"
    ),
    RENALFUNC_REL = list(
      description = "Relative renal function as a fraction of normal",
      units = "(dimensionless)",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "1 = normal renal function. Methods Section 4.2 maps KDIGO",
        "categories onto this scalar: normal = 1.0, mild impairment = 0.69,",
        "moderate = 0.32 and severe = 0.19. It multiplies ONLY the urinary",
        "excretion of M1 and M2 (unchanged glimepiride is not excreted",
        "renally), so renal impairment raises metabolite exposure and leaves",
        "parent exposure unchanged. Source name f_renal_function (SBML",
        "KI__f_renal_function). The authors' archived scenario simulations",
        "use 0.495, 0.347 and 0.198 for the three impairment classes",
        "instead; see the vignette Errata."
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
        "f_cirrhosis on [0, 0.95] (0 = normal), so this column is the",
        "complement: f_cirrhosis = 1 - HEPFUNC_REL. Methods Section 4.2",
        "maps f_cirrhosis onto the Child-Turcotte-Pugh classes following",
        "Koeller et al. (2021); the authors' archived scenario simulations",
        "use f_cirrhosis = 0.3995 (CTP A), 0.6980 (CTP B) and 0.8128",
        "(CTP C), i.e. HEPFUNC_REL = 0.6005, 0.3020 and 0.1872.",
        "f_cirrhosis acts twice in the source model, as the fraction of",
        "parenchymal liver volume lost (shrinking the metabolising tissue",
        "volume) and as the fraction of portal and hepatic-arterial blood",
        "shunted past the liver."
      ),
      source_name = "f_cirrhosis"
    ),
    CYP2C9 = list(
      description = "CYP2C9 metabolic activity relative to the *1/*1 genotype",
      units = "(fraction of *1/*1 activity)",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "1 = wild-type (*1/*1) activity. Multiplies the maximal rate of",
        "the CYP2C9-mediated glimepiride -> M1 conversion in liver tissue",
        "(source parameter f_cyp2c9, SBML LI__f_cyp2c9). Methods Section",
        "4.2 and Supplementary Table S2: allele activities *1 = 1.00,",
        "*2 = 0.68, *3 = 0.23; the diplotype value is the MEAN of its two",
        "allele activities, giving *1/*1 = 1.00, *1/*2 = 0.84,",
        "*1/*3 = 0.615 (printed 0.62), *2/*2 = 0.68, *2/*3 = 0.455 and",
        "*3/*3 = 0.23. For the paper's probabilistic population",
        "simulations each allele activity was additionally drawn from a",
        "lognormal distribution with shape sigma = 0.81 scaled to the",
        "allele mean (Supplementary Tables S4-S5)."
      ),
      source_name = "f_cyp2c9"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "glimepiride", units = "mg", specimen = "administration site", verified = TRUE),
    depot_iv = list(analyte = "glimepiride", units = "mg", specimen = "administration site", verified = TRUE),
    depot_iv_m1 = list(analyte = "M1", units = "mg", specimen = "administration site", verified = TRUE),
    depot_iv_m2 = list(analyte = "M2", units = "mg", specimen = "administration site", verified = TRUE),
    gut_lumen = list(analyte = "glimepiride", units = "mmol/L", specimen = "administration site", verified = TRUE),
    gut_lumen_m1 = list(analyte = "M1", units = "mmol/L", specimen = "faeces", verified = TRUE),
    gut_lumen_m2 = list(analyte = "M2", units = "mmol/L", specimen = "faeces", verified = TRUE),
    feces_m1 = list(analyte = "M1", units = "mmol", specimen = "faeces", verified = TRUE),
    feces_m2 = list(analyte = "M2", units = "mmol", specimen = "faeces", verified = TRUE),
    gut_plasma = list(analyte = "glimepiride", units = "mmol/L", specimen = "plasma", verified = TRUE),
    gut_plasma_m1 = list(analyte = "M1", units = "mmol/L", specimen = "plasma", verified = TRUE),
    gut_plasma_m2 = list(analyte = "M2", units = "mmol/L", specimen = "plasma", verified = TRUE),
    portal = list(analyte = "glimepiride", units = "mmol/L", specimen = "plasma", verified = TRUE),
    portal_m1 = list(analyte = "M1", units = "mmol/L", specimen = "plasma", verified = TRUE),
    portal_m2 = list(analyte = "M2", units = "mmol/L", specimen = "plasma", verified = TRUE),
    liver_plasma = list(analyte = "glimepiride", units = "mmol/L", specimen = "plasma", verified = TRUE),
    liver_plasma_m1 = list(analyte = "M1", units = "mmol/L", specimen = "plasma", verified = TRUE),
    liver_plasma_m2 = list(analyte = "M2", units = "mmol/L", specimen = "plasma", verified = TRUE),
    liver = list(analyte = "glimepiride", units = "mmol/L", specimen = "tissue", verified = TRUE),
    liver_m1 = list(analyte = "M1", units = "mmol/L", specimen = "tissue", verified = TRUE),
    liver_m2 = list(analyte = "M2", units = "mmol/L", specimen = "tissue", verified = TRUE),
    hepatic_vein = list(analyte = "glimepiride", units = "mmol/L", specimen = "plasma", verified = TRUE),
    hepatic_vein_m1 = list(analyte = "M1", units = "mmol/L", specimen = "plasma", verified = TRUE),
    hepatic_vein_m2 = list(analyte = "M2", units = "mmol/L", specimen = "plasma", verified = TRUE),
    kidney_plasma = list(analyte = "glimepiride", units = "mmol/L", specimen = "plasma", verified = TRUE),
    kidney_plasma_m1 = list(analyte = "M1", units = "mmol/L", specimen = "plasma", verified = TRUE),
    kidney_plasma_m2 = list(analyte = "M2", units = "mmol/L", specimen = "plasma", verified = TRUE),
    urine_m1 = list(analyte = "M1", units = "mmol", specimen = "urine", verified = TRUE),
    urine_m2 = list(analyte = "M2", units = "mmol", specimen = "urine", verified = TRUE),
    lung_plasma = list(analyte = "glimepiride", units = "mmol/L", specimen = "plasma", verified = TRUE),
    lung_plasma_m1 = list(analyte = "M1", units = "mmol/L", specimen = "plasma", verified = TRUE),
    lung_plasma_m2 = list(analyte = "M2", units = "mmol/L", specimen = "plasma", verified = TRUE),
    lung = list(analyte = "glimepiride", units = "mmol/L", specimen = "tissue", verified = TRUE),
    lung_m1 = list(analyte = "M1", units = "mmol/L", specimen = "tissue", verified = TRUE),
    lung_m2 = list(analyte = "M2", units = "mmol/L", specimen = "tissue", verified = TRUE),
    rest_plasma = list(analyte = "glimepiride", units = "mmol/L", specimen = "plasma", verified = TRUE),
    rest_plasma_m1 = list(analyte = "M1", units = "mmol/L", specimen = "plasma", verified = TRUE),
    rest_plasma_m2 = list(analyte = "M2", units = "mmol/L", specimen = "plasma", verified = TRUE),
    rest = list(analyte = "glimepiride", units = "mmol/L", specimen = "tissue", verified = TRUE),
    rest_m1 = list(analyte = "M1", units = "mmol/L", specimen = "tissue", verified = TRUE),
    rest_m2 = list(analyte = "M2", units = "mmol/L", specimen = "tissue", verified = TRUE),
    venous = list(analyte = "glimepiride", units = "mmol/L", specimen = "plasma", verified = TRUE),
    venous_m1 = list(analyte = "M1", units = "mmol/L", specimen = "plasma", verified = TRUE),
    venous_m2 = list(analyte = "M2", units = "mmol/L", specimen = "plasma", verified = TRUE),
    arterial = list(analyte = "glimepiride", units = "mmol/L", specimen = "plasma", verified = TRUE),
    arterial_m1 = list(analyte = "M1", units = "mmol/L", specimen = "plasma", verified = TRUE),
    arterial_m2 = list(analyte = "M2", units = "mmol/L", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = NA_integer_,
    n_studies = 20L,
    age_range = "adults; pediatric studies excluded (Methods Section 4.1)",
    weight_range = "75 kg reference individual; 40-170 kg explored in the body-weight scan (Figure 5)",
    sex_female_pct = NA_real_,
    disease_state = paste(
      "Healthy volunteers, patients with type 2 diabetes mellitus, subjects",
      "with renal impairment, normal-weight and morbidly obese subjects, and",
      "carriers of CYP2C9 *1, *2 and *3 alleles (Table 1). Hepatic",
      "impairment was simulated from the cirrhosis model with limited",
      "clinical comparison data."
    ),
    dose_range = paste(
      "Oral 0.5-8 mg single and multiple dose; intravenous glimepiride",
      "0.25-1.5 mg and intravenous M1 1.5 mg (Table 1)."
    ),
    regions = "Europe, East Asia, Egypt and the USA across the curated studies (Table 1 and Supplementary Section S6)",
    notes = paste(
      "Deterministic typical-individual model. The paper reports a single",
      "optimized parameter set (Supplementary Table S1), used for every",
      "simulation without refitting, and no inter-individual variability",
      "or residual-error model; the only stochastic element in the paper",
      "is the lognormal sampling of CYP2C9 activity in the population",
      "analyses. Parameters were fitted by weighted least squares (100",
      "local-optimisation runs) to the healthy, fasted subset of the",
      "curated PK-DB dataset (Methods Section 4.3)."
    )
  )

  ini({
    # ================= Reference physiology ==============================
    # The whole-body physiology is not printed in the paper; values are the
    # parameters of the authors' archived SBML v0.6.1 (cited in Methods 4.2).
    covbw   <- fixed(1.548);  label("Cardiac output per unit body weight (mL/s/kg)")  # SBML COBW
    fcardio <- fixed(1);      label("Cardiac function scaling factor (unitless)")     # SBML f_cardiac_function, normal = 1
    hct     <- fixed(0.51);   label("Hematocrit (fraction)")                          # SBML HCT
    fblood  <- fixed(0.02);   label("Blood fraction of organ volume (unitless)")      # SBML Fblood

    fvar <- fixed(0.0257);  label("Arterial blood fractional volume (L/kg)")    # SBML FVar
    fvgu <- fixed(0.0171);  label("Gut fractional volume (L/kg)")               # SBML FVgu
    fvhv <- fixed(0.001);   label("Hepatic venous fractional volume (L/kg)")    # SBML FVhv
    fvki <- fixed(0.0044);  label("Kidney fractional volume (L/kg)")            # SBML FVki
    fvli <- fixed(0.021);   label("Liver fractional volume (L/kg)")             # SBML FVli
    fvlu <- fixed(0.0076);  label("Lung fractional volume (L/kg)")              # SBML FVlu
    fvpo <- fixed(0.001);   label("Portal fractional volume (L/kg)")            # SBML FVpo
    fvve <- fixed(0.0514);  label("Venous blood fractional volume (L/kg)")      # SBML FVve

    fqgu <- fixed(0.18);   label("Gut fractional blood flow (unitless)")                  # SBML FQgu
    fqh  <- fixed(0.215);  label("Hepatic venous-side fractional blood flow (unitless)")  # SBML FQh
    fqki <- fixed(0.19);   label("Kidney fractional blood flow (unitless)")               # SBML FQki
    fqlu <- fixed(1);      label("Lung fractional blood flow (unitless)")                 # SBML FQlu

    # ================= Molecular weights =================================
    mr_gli <- fixed(490.616); label("Molecular weight of glimepiride (g/mol)")  # SBML Mr_gli (= GU__Mr_gli); Supplementary S3 'Mr gli'
    mr_m1  <- fixed(506.62);  label("Molecular weight of M1 (g/mol)")           # SBML Mr_m1
    mr_m2  <- fixed(520.6);   label("Molecular weight of M2 (g/mol)")           # SBML Mr_m2

    # ================= Distribution ======================================
    # Supplementary Table S1 (optimized). Methods 4.2: one partition
    # coefficient and one distribution rate shared by glimepiride, M1 and M2.
    kp_gli      <- fixed(10.0206);  label("Tissue-to-plasma partition coefficient, glimepiride, M1 and M2 (unitless)")  # Table S1 Kp gli = 10.02060
    ftissue_gli <- fixed(7.0864e-4); label("Tissue distribution rate, glimepiride, M1 and M2 (L/min)")                   # Table S1 ftissue gli = 0.00071 (printed to 2 significant digits; SBML 0.000708642)

    # ================= Absorption and intestinal metabolite handling =====
    ka_dis_gli   <- fixed(2);        label("Dissolution rate of glimepiride in the GI tract (1/hr)")             # SBML GU__Ka_dis_gli; Supplementary S3 'Ka dis_gli [hr-1]' (value not printed)
    gliabs_k     <- fixed(0.01590);  label("Absorption rate of glimepiride from the gut lumen (1/min)")           # Table S1 GU__GLIABS k = 0.01590
    f_absorption <- fixed(1);        label("Absorption-rate scaling factor for food effects (unitless)")         # SBML GU__f_absorption; Supplementary S3 'f absorption', fasted = 1
    mreabs_k     <- fixed(0.01592);  label("Rate of M1 and M2 transfer from gut plasma into the gut lumen (1/min)")  # Table S1 GU__MREABS k = 0.01592
    mexc_k       <- fixed(1.7193e-4); label("Fecal excretion rate of M1 and M2 from the gut lumen (1/min)")       # Table S1 GU__MEXC k = 0.00017 (printed to 2 significant digits; SBML 0.000171930)

    # ================= Hepatic processes =================================
    gliim_k        <- fixed(100);      label("Glimepiride import rate into liver tissue (1/min)")               # SBML LI__GLIIM_k (not optimized; Supplementary S3 'GLIIM k')
    gli2m1_vmax    <- fixed(5.362e-5); label("Vmax, CYP2C9 conversion of glimepiride to M1 (mmol/min/L)")      # Table S1 LI__GLI2M1 Vmax = 0.00005 (printed to 1 significant digit; SBML 5.36204e-5)
    gli2m1_km_gli  <- fixed(0.00225);  label("Km, CYP2C9 conversion of glimepiride to M1 (mmol/L)")            # SBML LI__GLI2M1_Km_gli; Methods 4.2 literature-parameterised (Suzuki 2006, Maekawa 2009, Zhang 2023)
    m1ex_li_k      <- fixed(0.07774);  label("M1 export rate from liver to plasma (1/min)")                    # Table S1 LI__M1EX k = 0.07774
    m12m2_k        <- fixed(0.01485);  label("M1 to M2 conversion rate in liver (1/min)")                      # Table S1 LI__M12M2 k = 0.01485
    m2ex_li_k      <- fixed(99.99817); label("M2 export rate from liver to plasma (1/min)")                    # Table S1 LI__M2EX k = 99.99817 (at its optimisation bound, Discussion)

    # ================= Renal excretion ===================================
    m1ex_ki_k <- fixed(0.14801); label("Urinary excretion rate of M1 (1/min)")  # Table S1 KI__M1EX k = 0.14801
    m2ex_ki_k <- fixed(0.09849); label("Urinary excretion rate of M2 (1/min)")  # Table S1 KI__M2EX k = 0.09849

    # ================= Intravenous input =================================
    ti_iv <- fixed(10); label("Nominal intravenous injection time (s)")  # SBML ti_gli = ti_m1 = ti_m2; drives a first-order release with half-life ti_iv

    # ================= Residual error ====================================
    # The paper reports a deterministic typical-individual model with no
    # residual-error model; fixed at zero rather than invented.
    propSd <- fixed(0); label("Proportional residual error (fraction; ZERO - not reported by the source)")  # Methods 4.3: single optimized parameter set, no error model
  })

  model({
    # ================= Scenario scalars ==================================
    # Renal function scales M1 and M2 urinary excretion only (Methods 4.2).
    f_renal <- RENALFUNC_REL

    # The source parameterises hepatic impairment as cirrhosis SEVERITY,
    # the complement of the canonical fraction-of-normal column.
    f_cirrhosis   <- 1 - HEPFUNC_REL
    f_shunts      <- f_cirrhosis
    f_tissue_loss <- f_cirrhosis

    # CYP2C9 activity relative to *1/*1 scales the M1-formation Vmax.
    f_cyp2c9 <- CYP2C9

    # ================= Volumes ===========================================
    fvre <- 1 - (fvgu + fvki + fvli + fvlu + fvve + fvar)

    vgu <- WT * fvgu
    vki <- WT * fvki
    vli <- WT * fvli
    vlu <- WT * fvlu
    vre <- WT * fvre

    # Arterial and venous blood pools are corrected for the blood held
    # within the organ volumes.
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

    # ================= Dosing compartments ===============================
    # Oral drug (mg) dissolves first-order into the gut lumen (mmol);
    # intravenous doses (mg) empty by a fast first-order release.
    dissolution <- (ka_dis_gli / 60) * depot / mr_gli

    ki_iv <- (0.693 / ti_iv) * 60
    iv_gli <- ki_iv * depot_iv / mr_gli
    iv_m1  <- ki_iv * depot_iv_m1 / mr_m1
    iv_m2  <- ki_iv * depot_iv_m2 / mr_m2

    # ================= Intestine (Supplementary S3) ======================
    gliabs <- f_absorption * gliabs_k * vgu * gut_lumen
    m1reabs <- mreabs_k * gut_plasma_m1 * vgu
    m2reabs <- mreabs_k * gut_plasma_m2 * vgu
    m1exc <- mexc_k * gut_lumen_m1 * vgu
    m2exc <- mexc_k * gut_lumen_m2 * vgu

    # ================= Liver (Supplementary S3) ==========================
    li_gliim  <- gliim_k * vli_tissue * (liver_plasma - liver)
    li_gli2m1 <- f_cyp2c9 * gli2m1_vmax * vli_tissue * liver / (liver + gli2m1_km_gli)
    li_m12m2  <- m12m2_k * vli_tissue * liver_m1
    li_m1ex   <- m1ex_li_k * vli_tissue * (liver_m1 - liver_plasma_m1)
    li_m2ex   <- m2ex_li_k * vli_tissue * (liver_m2 - liver_plasma_m2)

    # ================= Kidney (Supplementary S3) =========================
    # Excretion is first order in kidney PLASMA concentration, scaled by the
    # kidney tissue volume, as in the source.
    ki_m1ex <- f_renal * vki_tissue * m1ex_ki_k * kidney_plasma_m1
    ki_m2ex <- f_renal * vki_tissue * m2ex_ki_k * kidney_plasma_m2

    # ================= Tissue distribution (lung, rest of body) ==========
    transport_lu    <- ftissue_gli * (lung_plasma * kp_gli - lung)
    transport_lu_m1 <- ftissue_gli * (lung_plasma_m1 * kp_gli - lung_m1)
    transport_lu_m2 <- ftissue_gli * (lung_plasma_m2 * kp_gli - lung_m2)
    transport_re    <- ftissue_gli * (rest_plasma * kp_gli - rest)
    transport_re_m1 <- ftissue_gli * (rest_plasma_m1 * kp_gli - rest_m1)
    transport_re_m2 <- ftissue_gli * (rest_plasma_m2 * kp_gli - rest_m2)

    # ================= ODEs: dosing ======================================
    d/dt(depot)       <- -dissolution * mr_gli
    d/dt(depot_iv)    <- -iv_gli * mr_gli
    d/dt(depot_iv_m1) <- -iv_m1 * mr_m1
    d/dt(depot_iv_m2) <- -iv_m2 * mr_m2

    # ================= ODEs: gut lumen and feces =========================
    d/dt(gut_lumen)    <- -gliabs / vgu + dissolution / vgu
    d/dt(gut_lumen_m1) <- m1reabs / vgu - m1exc / vgu
    d/dt(gut_lumen_m2) <- m2reabs / vgu - m2exc / vgu
    d/dt(feces_m1)     <- m1exc
    d/dt(feces_m2)     <- m2exc

    # ================= ODEs: gut plasma ==================================
    d/dt(gut_plasma)    <- qgu * arterial / vgu_plasma - qgu * gut_plasma / vgu_plasma + gliabs / vgu_plasma
    d/dt(gut_plasma_m1) <- qgu * arterial_m1 / vgu_plasma - qgu * gut_plasma_m1 / vgu_plasma - m1reabs / vgu_plasma
    d/dt(gut_plasma_m2) <- qgu * arterial_m2 / vgu_plasma - qgu * gut_plasma_m2 / vgu_plasma - m2reabs / vgu_plasma

    # ================= ODEs: portal vein =================================
    d/dt(portal)    <- (qgu * gut_plasma - qpo * portal) / vpo
    d/dt(portal_m1) <- (qgu * gut_plasma_m1 - qpo * portal_m1) / vpo
    d/dt(portal_m2) <- (qgu * gut_plasma_m2 - qpo * portal_m2) / vpo

    # ================= ODEs: liver plasma and tissue =====================
    # Cirrhotic shunting diverts portal and hepatic-arterial blood straight
    # into the hepatic vein, bypassing the metabolising liver.
    d/dt(liver_plasma) <- ((1 - f_shunts) * qha * arterial + (1 - f_shunts) * qpo * portal -
      (1 - f_shunts) * (qpo + qha) * liver_plasma - li_gliim) / vli_plasma
    d/dt(liver_plasma_m1) <- ((1 - f_shunts) * qha * arterial_m1 + (1 - f_shunts) * qpo * portal_m1 -
      (1 - f_shunts) * (qpo + qha) * liver_plasma_m1 + li_m1ex) / vli_plasma
    d/dt(liver_plasma_m2) <- ((1 - f_shunts) * qha * arterial_m2 + (1 - f_shunts) * qpo * portal_m2 -
      (1 - f_shunts) * (qpo + qha) * liver_plasma_m2 + li_m2ex) / vli_plasma

    d/dt(liver)    <- (li_gliim - li_gli2m1) / vli_tissue
    d/dt(liver_m1) <- (li_gli2m1 - li_m1ex - li_m12m2) / vli_tissue
    d/dt(liver_m2) <- (li_m12m2 - li_m2ex) / vli_tissue

    # ================= ODEs: hepatic vein ================================
    d/dt(hepatic_vein) <- (f_shunts * qha * arterial + f_shunts * qpo * portal +
      (1 - f_shunts) * (qpo + qha) * liver_plasma - qh * hepatic_vein) / vhv
    d/dt(hepatic_vein_m1) <- (f_shunts * qha * arterial_m1 + f_shunts * qpo * portal_m1 +
      (1 - f_shunts) * (qpo + qha) * liver_plasma_m1 - qh * hepatic_vein_m1) / vhv
    d/dt(hepatic_vein_m2) <- (f_shunts * qha * arterial_m2 + f_shunts * qpo * portal_m2 +
      (1 - f_shunts) * (qpo + qha) * liver_plasma_m2 - qh * hepatic_vein_m2) / vhv

    # ================= ODEs: kidney and urine ============================
    d/dt(kidney_plasma)    <- (qki * arterial - qki * kidney_plasma) / vki_plasma
    d/dt(kidney_plasma_m1) <- (qki * arterial_m1 - qki * kidney_plasma_m1 - ki_m1ex) / vki_plasma
    d/dt(kidney_plasma_m2) <- (qki * arterial_m2 - qki * kidney_plasma_m2 - ki_m2ex) / vki_plasma
    d/dt(urine_m1) <- ki_m1ex
    d/dt(urine_m2) <- ki_m2ex

    # ================= ODEs: lung ========================================
    d/dt(lung_plasma)    <- (-transport_lu + qlu * venous - qlu * lung_plasma) / vlu_plasma
    d/dt(lung_plasma_m1) <- (-transport_lu_m1 + qlu * venous_m1 - qlu * lung_plasma_m1) / vlu_plasma
    d/dt(lung_plasma_m2) <- (-transport_lu_m2 + qlu * venous_m2 - qlu * lung_plasma_m2) / vlu_plasma
    d/dt(lung)    <- transport_lu / vlu_tissue
    d/dt(lung_m1) <- transport_lu_m1 / vlu_tissue
    d/dt(lung_m2) <- transport_lu_m2 / vlu_tissue

    # ================= ODEs: rest of body ================================
    d/dt(rest_plasma)    <- (-transport_re + qre * arterial - qre * rest_plasma) / vre_plasma
    d/dt(rest_plasma_m1) <- (-transport_re_m1 + qre * arterial_m1 - qre * rest_plasma_m1) / vre_plasma
    d/dt(rest_plasma_m2) <- (-transport_re_m2 + qre * arterial_m2 - qre * rest_plasma_m2) / vre_plasma
    d/dt(rest)    <- transport_re / vre_tissue
    d/dt(rest_m1) <- transport_re_m1 / vre_tissue
    d/dt(rest_m2) <- transport_re_m2 / vre_tissue

    # ================= ODEs: systemic circulation ========================
    d/dt(venous) <- (iv_gli + qki * kidney_plasma + qh * hepatic_vein -
      qlu * venous + qre * rest_plasma) / vve
    d/dt(venous_m1) <- (iv_m1 + qki * kidney_plasma_m1 + qh * hepatic_vein_m1 -
      qlu * venous_m1 + qre * rest_plasma_m1) / vve
    d/dt(venous_m2) <- (iv_m2 + qki * kidney_plasma_m2 + qh * hepatic_vein_m2 -
      qlu * venous_m2 + qre * rest_plasma_m2) / vve

    d/dt(arterial) <- (-qki * arterial - (1 - f_shunts) * qha * arterial - f_shunts * qha * arterial +
      qlu * lung_plasma - qgu * arterial - qre * arterial) / var_
    d/dt(arterial_m1) <- (-qki * arterial_m1 - (1 - f_shunts) * qha * arterial_m1 - f_shunts * qha * arterial_m1 +
      qlu * lung_plasma_m1 - qgu * arterial_m1 - qre * arterial_m1) / var_
    d/dt(arterial_m2) <- (-qki * arterial_m2 - (1 - f_shunts) * qha * arterial_m2 - f_shunts * qha * arterial_m2 +
      qlu * lung_plasma_m2 - qgu * arterial_m2 - qre * arterial_m2) / var_

    # ================= Observations ======================================
    # Venous plasma concentrations in umol/L, matching the paper's figures.
    Cc    <- 1000 * venous
    Cc_m1 <- 1000 * venous_m1
    Cc_m2 <- 1000 * venous_m2

    # Cumulative urinary and fecal recovery of the metabolites (mmol).
    Aurine_m1m2 <- urine_m1 + urine_m2
    Afeces_m1m2 <- feces_m1 + feces_m2

    Cc ~ prop(propSd)
  })
}
