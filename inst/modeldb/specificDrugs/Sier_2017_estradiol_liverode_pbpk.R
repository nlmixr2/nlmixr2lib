Sier_2017_estradiol_liverode_pbpk <- function() {
  description <- paste(
    "PBPK (whole-body, COPASI). Model of endogenous and exogenous",
    "17-beta-estradiol (E2)",
    "in women with a mechanistic liver-metabolism network in place of a",
    "single intrinsic-clearance term (PBPK-LiverODE variant of Sier 2017).",
    "Same skeleton as modellib('Sier_2017_estradiol_pbpk') -- venous and",
    "arterial blood with SHBG and albumin binding, 16 tissues,",
    "permeability-limited liver and uterus+ovaries 'gonads', estrogen-receptor",
    "binding, zero-order gonadal biosynthesis, kidney extra-hepatic clearance",
    "-- but hepatocyte E2 is cleared by 27 Michaelis-Menten enzymatic",
    "conversions (UGT1A1/1A3/2B7, SULT1E1/2A1, CYP1A2/3A4, COMT, HSD17B1/2)",
    "through estrone and catechol intermediates, with non-specific binding and",
    "6 efflux transporter reactions (ABCG2, MRP2/3/4, MDR1) for the conjugates.",
    "Enzyme amounts are MOPED protein abundances. Blood flows are the",
    "deposited volume * specific flow * 0.55 (no tissue density), with lung",
    "flow 0.55 * cardiac output, and liver permeability 150 / 41.7 1/h, all as",
    "deposited for this variant. States start at the deposited endogenous",
    "steady state of the 60 kg, 1.63 m, 35-year reference woman. Observation",
    "Cc is FREE (unbound) venous E2."
  )
  reference <- paste(
    "Sier JH, Thumser AE, Plant NJ. Linking physiologically-based",
    "pharmacokinetic and genome-scale metabolic networks to understand",
    "estradiol biology. BMC Syst Biol. 2017;11:141.",
    "doi:10.1186/s12918-017-0520-3. Executable model deposited by the authors",
    "as BioModels MODEL1711210002 (E2_PBPK_LiverODE_v7.3)."
  )
  vignette <- "Sier_2017_estradiol"
  # Every ODE state is a concentration (nM); the deposited COPASI model sets
  # every compartment size to 1, so a dose in nmol entered into a dose
  # compartment is numerically its starting concentration.
  units <- list(time = "h", dosing = "nmol", concentration = "nM")

  dosing <- c("depot", "depot_iv")

  # Hepatocyte species of the metabolism network. They are specific to this
  # estradiol network (estrone, catechol estrogens, conjugates and a
  # non-specific binding pool), so they are whitelisted here rather than
  # registered as canonical compartments.
  paper_specific_compartments <- c(
    "e1_liver",
    "oh2e2_liver",
    "oh4e2_liver",
    "oh2e1_liver",
    "oh4e1_liver",
    "e13s_liver",
    "e217bg_liver",
    "nsb_liver",
    "complex_nsb_liver"
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Drives body surface area, body mass index, adipose volume and,",
        "through the 93% lean-mass balance, muscle volume. Deposited",
        "reference 60 kg."
      ),
      source_name = "Body weight (kg)"
    ),
    HT = list(
      description = "Body height",
      units = "cm",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "The source equations take height in metres; model() converts",
        "ht_m = HT / 100. Deposited reference 1.63 m (163 cm)."
      ),
      source_name = "Height (m)"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Drives adipose, brain, gonad and thymus volumes. Deposited reference 35 years.",
      source_name = "Age (years)"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "estradiol", units = "nmol", specimen = "administration site", verified = TRUE),
    depot_iv = list(analyte = "estradiol", units = "nmol", specimen = "administration site", verified = TRUE),
    venous = list(analyte = "estradiol (free)", units = "nM", specimen = "plasma", verified = TRUE),
    shbg_venous = list(analyte = "SHBG (free)", units = "nM", specimen = "plasma", verified = TRUE),
    complex_shbg_venous = list(analyte = "estradiol-SHBG complex", units = "nM", specimen = "plasma", verified = TRUE),
    albumin_venous = list(analyte = "albumin (free)", units = "nM", specimen = "plasma", verified = TRUE),
    complex_albumin_venous = list(
      analyte = "estradiol-albumin complex",
      units = "nM",
      specimen = "plasma",
      verified = TRUE
    ),
    arterial = list(analyte = "estradiol (free)", units = "nM", specimen = "plasma", verified = TRUE),
    shbg_arterial = list(analyte = "SHBG (free)", units = "nM", specimen = "plasma", verified = TRUE),
    complex_shbg_arterial = list(
      analyte = "estradiol-SHBG complex",
      units = "nM",
      specimen = "plasma",
      verified = TRUE
    ),
    albumin_arterial = list(analyte = "albumin (free)", units = "nM", specimen = "plasma", verified = TRUE),
    complex_albumin_arterial = list(
      analyte = "estradiol-albumin complex",
      units = "nM",
      specimen = "plasma",
      verified = TRUE
    ),
    lung = list(analyte = "estradiol", units = "nM", specimen = "tissue", verified = TRUE),
    adipose = list(analyte = "estradiol", units = "nM", specimen = "tissue", verified = TRUE),
    bone = list(analyte = "estradiol", units = "nM", specimen = "tissue", verified = TRUE),
    brain = list(analyte = "estradiol", units = "nM", specimen = "tissue", verified = TRUE),
    heart = list(analyte = "estradiol", units = "nM", specimen = "tissue", verified = TRUE),
    muscle = list(analyte = "estradiol", units = "nM", specimen = "tissue", verified = TRUE),
    skin = list(analyte = "estradiol", units = "nM", specimen = "tissue", verified = TRUE),
    intestine = list(analyte = "estradiol", units = "nM", specimen = "tissue", verified = TRUE),
    spleen = list(analyte = "estradiol", units = "nM", specimen = "tissue", verified = TRUE),
    stomach = list(analyte = "estradiol", units = "nM", specimen = "tissue", verified = TRUE),
    pancreas = list(analyte = "estradiol", units = "nM", specimen = "tissue", verified = TRUE),
    thymus = list(analyte = "estradiol", units = "nM", specimen = "tissue", verified = TRUE),
    kidney = list(analyte = "estradiol", units = "nM", specimen = "tissue", verified = TRUE),
    other = list(analyte = "estradiol", units = "nM", specimen = "tissue", verified = TRUE),
    vp_liver = list(analyte = "estradiol", units = "nM", specimen = "whole blood", verified = TRUE),
    liver = list(analyte = "estradiol (free)", units = "nM", specimen = "tissue", verified = TRUE),
    target_er_liver = list(analyte = "estrogen receptor (free)", units = "nM", specimen = "tissue", verified = TRUE),
    complex_er_liver = list(
      analyte = "estradiol-estrogen receptor complex",
      units = "nM",
      specimen = "tissue",
      verified = TRUE
    ),
    nsb_liver = list(analyte = "non-specific binding sites (free)", units = "nM", specimen = "tissue", verified = TRUE),
    complex_nsb_liver = list(
      analyte = "estradiol bound to non-specific sites",
      units = "nM",
      specimen = "tissue",
      verified = TRUE
    ),
    e1_liver = list(analyte = "estrone", units = "nM", specimen = "tissue", verified = TRUE),
    oh2e2_liver = list(analyte = "2-hydroxyestradiol", units = "nM", specimen = "tissue", verified = TRUE),
    oh4e2_liver = list(analyte = "4-hydroxyestradiol", units = "nM", specimen = "tissue", verified = TRUE),
    oh2e1_liver = list(analyte = "2-hydroxyestrone", units = "nM", specimen = "tissue", verified = TRUE),
    oh4e1_liver = list(analyte = "4-hydroxyestrone", units = "nM", specimen = "tissue", verified = TRUE),
    e13s_liver = list(analyte = "estrone-3-sulfate", units = "nM", specimen = "tissue", verified = TRUE),
    e217bg_liver = list(analyte = "estradiol-17beta-glucuronide", units = "nM", specimen = "tissue", verified = TRUE),
    vp_gonads = list(analyte = "estradiol", units = "nM", specimen = "whole blood", verified = TRUE),
    gonads = list(analyte = "estradiol (free)", units = "nM", specimen = "tissue", verified = TRUE),
    target_er_gonads = list(analyte = "estrogen receptor (free)", units = "nM", specimen = "tissue", verified = TRUE),
    complex_er_gonads = list(
      analyte = "estradiol-estrogen receptor complex",
      units = "nM",
      specimen = "tissue",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 1495L,
    n_studies = 2L,
    age_range = "18-45 years (virtual population)",
    sex_female_pct = 100,
    disease_state = "healthy women (endogenous estradiol plus single exogenous doses)",
    dose_range = "Single 0.3 mg i.v. bolus; single 2 mg or 4 mg oral dose",
    regions = "United States (NHANES anthropometrics)",
    notes = paste(
      "Not a population fit: a deterministic mechanistic model with no",
      "inter-individual variability or residual error. The 1495 subjects are",
      "a virtual population of women aged 18-45 years drawn from NHANES",
      "2013-14 anthropometrics. Compared against Kuhnz 1993 (0.3 mg i.v.)",
      "and Lyrenas 1981 (2 and 4 mg oral). The stand-alone liver network was",
      "checked against estradiol clearance by primary human hepatocytes",
      "(Additional file 1 Figure S4)."
    )
  )

  ini({
    # All values are point values from the paper / Additional file 1 / the
    # authors' BioModels deposit MODEL1711210002; nothing was estimated.

    # E2 input -----------------------------------------------------------
    ksyn_gonads <- fixed(2)
    label("E2 biosynthesis rate into gonad cells (nM/h)") # BioModels MODEL1711210002 'Gonads: E2 Biosynthesis' k1 = 2 times the 1 nM Biosynthesis species
    kiv <- fixed(250)
    label("First-order delivery rate constant from the i.v. dose compartment to venous blood (1/h)") # Main text Results: 250 h-1; BioModels 'Dose_IV bolus' k1 = 250
    lka <- fixed(log(0.025))
    label("First-order absorption rate constant from the oral dose compartment to the intestine (1/h)") # Main text Results: 0.025 h-1; BioModels 'Dose_oral' k1 = 0.025
    kfeces <- fixed(0.01)
    label("First-order loss of unabsorbed oral dose (faecal excretion) (1/h)") # Main text Results: 0.01 h-1; BioModels 'Oral elimination' k1 = 0.01

    # Extra-hepatic clearance ---------------------------------------------------
    kcleh <- fixed(5)
    label("Extra-hepatic clearance rate constant acting on kidney E2 (1/h)") # Additional file 1 v10 k1 = 5 h-1; BioModels 'CLeh' k1 = 5

    # Permeability-limited tissues -------------------------------------------
    kperm_in_liver <- fixed(150)
    label("Liver uptake rate constant, vascular to hepatocyte (1/h)") # Main text Results: for PBPK-LiverODE, uptake 150 L.h-1; BioModels 'Liver: Transport, E2' k1 = 150
    kperm_out_liver <- fixed(41.7)
    label("Liver efflux rate constant, hepatocyte to vascular (1/h)") # Main text Results: for PBPK-LiverODE, efflux 41.7 L.h-1; BioModels k2 = 41.7
    kperm_in_gonads <- fixed(1000)
    label("Gonad uptake rate constant, vascular to cells (1/h)") # BioModels 'Gonads: Permeability-limited access' k1 = 1000
    kperm_out_gonads <- fixed(100)
    label("Gonad efflux rate constant, cells to vascular (1/h)") # BioModels 'Gonads: Permeability-limited access' k2 = 100

    # Plasma protein binding ----------------------------------------------------
    kon_shbg <- fixed(3600)
    label("E2-SHBG association rate constant (1/(nM*h))") # BioModels re2 / re4 Kon = 3600
    koff_shbg <- fixed(5400)
    label("E2-SHBG dissociation rate constant (1/h)") # BioModels re2 / re4 koff = 5400
    kon_albumin <- fixed(3600)
    label("E2-albumin association rate constant (1/(nM*h))") # BioModels re3 / re5 Kon = 3600
    koff_albumin <- fixed(6.12e7)
    label("E2-albumin dissociation rate constant (1/h)") # BioModels re3 / re5 koff = 61200000
    bl_shbg <- fixed(40)
    label("Total (free + E2-bound) SHBG concentration in venous and arterial plasma (nM)") # BioModels MODEL1711210002 initial concentrations: free SHBG + E2.SHBG = 40 nM
    bl_albumin <- fixed(500000)
    label("Total (free + E2-bound) albumin concentration in venous and arterial plasma (nM)") # BioModels MODEL1711210002 initial concentrations: free albumin + E2.albumin = 500000 nM

    # Receptor and non-specific binding ----------------------------------------
    kon_er_liver <- fixed(3600)
    label("E2-estrogen receptor association rate constant in liver cells (1/(nM*h))") # Additional file 1 LiverODE v31 k1 = 3600 h-1; BioModels 'Liver: Binding, ER binding' Kon = 3600
    koff_er_liver <- fixed(900)
    label("E2-estrogen receptor dissociation rate constant in liver cells (1/h)") # Additional file 1 LiverODE v31 k2 = 900 h-1; BioModels koff = 900
    kon_er_gonads <- fixed(3600)
    label("E2-estrogen receptor association rate constant in gonad cells (1/(nM*h))") # BioModels 'Gonads: ER binding' k1 = 3600
    koff_er_gonads <- fixed(900)
    label("E2-estrogen receptor dissociation rate constant in gonad cells (1/h)") # BioModels 'Gonads: ER binding' k2 = 900
    bl_er_liver <- fixed(46)
    label("Total (free + E2-bound) estrogen receptor in liver cells (nM)") # BioModels MODEL1711210002 initial concentrations: ER + ER_E2 (liver cells) = 45.99919 + 0.00081 = 46 nM
    bl_er_gonads <- fixed(140)
    label("Total (free + E2-bound) estrogen receptor in gonad cells (nM)") # BioModels MODEL1711210002 initial concentrations: ER + ER_E2 (gonad cells) = 0.76207 + 139.23793 = 140 nM
    kon_nsb <- fixed(3600)
    label("E2 association rate constant with hepatocyte non-specific binding sites (1/(nM*h))") # BioModels 'Liver: Binding, NSB binding' k1 = 3600 (Additional file 1 LiverODE v32 prints 0.1; see vignette)
    koff_nsb <- fixed(85500)
    label("E2 dissociation rate constant from hepatocyte non-specific binding sites (1/h)") # BioModels 'Liver: Binding, NSB binding' k2 = 85500 (Additional file 1 LiverODE v32 prints 0.1; see vignette)
    bl_nsb_liver <- fixed(13.5)
    label("Total (free + E2-bound) hepatocyte non-specific binding sites (nM)") # BioModels MODEL1711210002 initial concentrations: NSB + E2_NSB = 13.4999975 + 0.0000025 = 13.5 nM

    # Enzyme and transporter abundances -----------------------------------------
    f_protein_expression <- fixed(1)
    label("Protein expression factor scaling every enzyme abundance except ABCG2 (unitless)") # BioModels MODEL1711210002 Protein_Expression_Factor = 1 (the in vitro hepatocyte model of Additional file 1 uses 0.1)
    abund_comt <- fixed(596)
    label("COMT abundance in hepatocytes (nM)") # Additional file 1 LiverODE Initial Conditions: COMT (heps) 596*Protein expression factor nM
    abund_cyp1a2 <- fixed(13878)
    label("CYP1A2 abundance in hepatocytes (nM)") # Additional file 1 LiverODE Initial Conditions: CYP1A2 (heps) 13878*Protein expression factor nM
    abund_cyp3a4 <- fixed(10219)
    label("CYP3A4 abundance in hepatocytes (nM)") # Additional file 1 LiverODE Initial Conditions: CYP3A4 (heps) 10219*Protein expression factor nM
    abund_hsd17b1 <- fixed(541)
    label("HSD17B1 abundance in hepatocytes (nM)") # Additional file 1 LiverODE Initial Conditions: HSD17b1 (heps) 541*Protein expression factor nM
    abund_hsd17b2 <- fixed(1435)
    label("HSD17B2 abundance in hepatocytes (nM)") # Additional file 1 LiverODE Initial Conditions: HSD17b2 (heps) 1435*Protein expression factor nM
    abund_sult1e1 <- fixed(2932)
    label("SULT1E1 abundance in hepatocytes (nM)") # Additional file 1 LiverODE Initial Conditions: SULT1E1 (heps) 2932*Protein expression factor nM
    abund_sult2a1 <- fixed(64067)
    label("SULT2A1 abundance in hepatocytes (nM)") # Additional file 1 LiverODE Initial Conditions: SULT2A1 (heps) 64067*Protein expression factor nM
    abund_ugt1a1 <- fixed(1980)
    label("UGT1A1 abundance in hepatocytes (nM)") # Additional file 1 LiverODE Initial Conditions: UGT1A1 (heps) 1980*Protein expression factor nM
    abund_ugt1a3 <- fixed(2684)
    label("UGT1A3 abundance in hepatocytes (nM)") # Additional file 1 LiverODE Initial Conditions: UGT1A3 (heps) 2684*Protein expression factor nM
    abund_ugt2b7 <- fixed(17959)
    label("UGT2B7 abundance in hepatocytes (nM)") # Additional file 1 LiverODE Initial Conditions: UGT2B7 (heps) 17959*Protein expression factor nM
    abund_abcg2 <- fixed(83)
    label("ABCG2 abundance in hepatocytes (nM)") # Additional file 1 LiverODE Initial Conditions: ABCG2 (heps) 83*Protein expression factor nM; the deposit holds it at 83 without the factor
    abund_mdr1 <- fixed(223)
    label("MDR1 abundance in hepatocytes (nM)") # Additional file 1 LiverODE Initial Conditions: MDR1 (heps) 223*Protein expression factor nM
    abund_mrp2 <- fixed(153)
    label("MRP2 abundance in hepatocytes (nM)") # Additional file 1 LiverODE Initial Conditions: MRP2 (heps) 153*Protein expression factor nM
    abund_mrp3 <- fixed(362)
    label("MRP3 abundance in hepatocytes (nM)") # Additional file 1 LiverODE Initial Conditions: MRP3 (heps) 362*Protein expression factor nM
    abund_mrp4 <- fixed(100)
    label("MRP4 abundance in hepatocytes (nM)") # Additional file 1 LiverODE Initial Conditions: MRP4 (heps) 100*Protein expression factor nM

    # Michaelis-Menten kinetics: rate = S * enzyme * kcat / (S + Km), with
    # every value as deposited in MODEL1711210002 (Additional file 1 LiverODE
    # Reactions table gives the same value except where noted).
    kcat_ugt1a1_e2_e23g <- fixed(470.6)
    label("Turnover number, E2 -> E2-3-glucuronide by UGT1A1 (1/h)") # Additional file 1 v27 kcat; deposit Kcat = 470.6
    km_ugt1a1_e2_e23g <- fixed(23000)
    label("Michaelis constant, E2 -> E2-3-glucuronide by UGT1A1 (nM)") # Additional file 1 v27 Km; deposit km = 23000
    kcat_ugt1a3_e2_e23g <- fixed(145.5)
    label("Turnover number, E2 -> E2-3-glucuronide by UGT1A3 (1/h)") # Additional file 1 v28 kcat; deposit Kcat = 145.5
    km_ugt1a3_e2_e23g <- fixed(47000)
    label("Michaelis constant, E2 -> E2-3-glucuronide by UGT1A3 (nM)") # Additional file 1 v28 Km; deposit km = 47000
    kcat_ugt1a1_e2_e217bg <- fixed(212.5)
    label("Turnover number, E2 -> E2-17beta-glucuronide by UGT1A1 (1/h)") # Additional file 1 v24 kcat; deposit Kcat = 212.5
    km_ugt1a1_e2_e217bg <- fixed(35000)
    label("Michaelis constant, E2 -> E2-17beta-glucuronide by UGT1A1 (nM)") # Additional file 1 v24 Km; deposit km = 35000
    kcat_ugt1a3_e2_e217bg <- fixed(48.5)
    label("Turnover number, E2 -> E2-17beta-glucuronide by UGT1A3 (1/h)") # Additional file 1 v25 kcat; deposit Kcat = 48.5
    km_ugt1a3_e2_e217bg <- fixed(35000)
    label("Michaelis constant, E2 -> E2-17beta-glucuronide by UGT1A3 (nM)") # Additional file 1 v25 Km; deposit km = 35000
    kcat_ugt2b7_e2_e217bg <- fixed(102.7)
    label("Turnover number, E2 -> E2-17beta-glucuronide by UGT2B7 (1/h)") # Additional file 1 v26 kcat; deposit Kcat = 102.7
    km_ugt2b7_e2_e217bg <- fixed(11000)
    label("Michaelis constant, E2 -> E2-17beta-glucuronide by UGT2B7 (nM)") # Additional file 1 v26 Km; deposit km = 11000
    kcat_sult2a1_e2_e23s <- fixed(93.7)
    label("Turnover number, E2 -> E2-3-sulfate by SULT2A1 (1/h)") # Additional file 1 v29 kcat; deposit Kcat = 93.7
    km_sult2a1_e2_e23s <- fixed(14)
    label("Michaelis constant, E2 -> E2-3-sulfate by SULT2A1 (nM)") # Additional file 1 v29 Km; deposit km = 14
    kcat_cyp1a2_e2_oh2e2 <- fixed(660)
    label("Turnover number, E2 -> 2-hydroxy-E2 by CYP1A2 (1/h)") # Additional file 1 v18 kcat; deposit Kcat = 660
    km_cyp1a2_e2_oh2e2 <- fixed(20000)
    label("Michaelis constant, E2 -> 2-hydroxy-E2 by CYP1A2 (nM)") # Additional file 1 v18 Km; deposit km = 20000
    kcat_cyp3a4_e2_oh2e2 <- fixed(48)
    label("Turnover number, E2 -> 2-hydroxy-E2 by CYP3A4 (1/h)") # Additional file 1 v19 kcat; deposit Kcat = 48
    km_cyp3a4_e2_oh2e2 <- fixed(54000)
    label("Michaelis constant, E2 -> 2-hydroxy-E2 by CYP3A4 (nM)") # Additional file 1 v19 Km; deposit km = 54000
    kcat_cyp1a2_e2_oh4e2 <- fixed(54)
    label("Turnover number, E2 -> 4-hydroxy-E2 by CYP1A2 (1/h)") # Additional file 1 v20 kcat; deposit Kcat = 54
    km_cyp1a2_e2_oh4e2 <- fixed(28000)
    label("Michaelis constant, E2 -> 4-hydroxy-E2 by CYP1A2 (nM)") # Additional file 1 v20 Km; deposit km = 28000
    kcat_cyp3a4_e2_oh4e2 <- fixed(18)
    label("Turnover number, E2 -> 4-hydroxy-E2 by CYP3A4 (1/h)") # Additional file 1 v21 kcat; deposit Kcat = 18
    km_cyp3a4_e2_oh4e2 <- fixed(111000)
    label("Michaelis constant, E2 -> 4-hydroxy-E2 by CYP3A4 (nM)") # Additional file 1 v21 Km; deposit km = 111000
    kcat_cyp1a2_e2_oh16ae2 <- fixed(42)
    label("Turnover number, E2 -> 16alpha-hydroxy-E2 by CYP1A2 (1/h)") # Additional file 1 v16 kcat; deposit Kcat = 42
    km_cyp1a2_e2_oh16ae2 <- fixed(58000)
    label("Michaelis constant, E2 -> 16alpha-hydroxy-E2 by CYP1A2 (nM)") # Additional file 1 v16 Km; deposit km = 58000
    kcat_cyp3a4_e2_oh16ae2 <- fixed(24)
    label("Turnover number, E2 -> 16alpha-hydroxy-E2 by CYP3A4 (1/h)") # Additional file 1 v17 kcat; deposit Kcat = 24
    km_cyp3a4_e2_oh16ae2 <- fixed(75000)
    label("Michaelis constant, E2 -> 16alpha-hydroxy-E2 by CYP3A4 (nM)") # Additional file 1 v17 Km; deposit km = 75000
    kcat_hsd17b1_e2_e1 <- fixed(43.4)
    label("Turnover number, E2 -> E1 by HSD17B1 (1/h)") # Additional file 1 v22 kcat; deposit Kcat = 43.4
    km_hsd17b1_e2_e1 <- fixed(2000)
    label("Michaelis constant, E2 -> E1 by HSD17B1 (nM)") # Additional file 1 v22 Km; deposit km = 2000
    kcat_hsd17b2_e2_e1 <- fixed(16926)
    label("Turnover number, E2 -> E1 by HSD17B2 (1/h)") # Additional file 1 v23 kcat; deposit Kcat = 16926
    km_hsd17b2_e2_e1 <- fixed(110000)
    label("Michaelis constant, E2 -> E1 by HSD17B2 (nM)") # Additional file 1 v23 Km; deposit km = 110000
    kcat_hsd17b1_e1_e2 <- fixed(5400)
    label("Turnover number, E1 -> E2 by HSD17B1 (1/h)") # Additional file 1 v15 kcat; deposit Kcat = 5400
    km_hsd17b1_e1_e2 <- fixed(70000)
    label("Michaelis constant, E1 -> E2 by HSD17B1 (nM)") # Additional file 1 v15 Km; deposit km = 70000
    kcat_ugt1a1_e1_e13g <- fixed(15.2)
    label("Turnover number, E1 -> E1-3-glucuronide by UGT1A1 (1/h)") # Additional file 1 v12 kcat; deposit Kcat = 15.2
    km_ugt1a1_e1_e13g <- fixed(38000)
    label("Michaelis constant, E1 -> E1-3-glucuronide by UGT1A1 (nM)") # Additional file 1 v12 Km; deposit km = 38000
    kcat_ugt1a3_e1_e13g <- fixed(186.6)
    label("Turnover number, E1 -> E1-3-glucuronide by UGT1A3 (1/h)") # Additional file 1 v13 kcat; deposit Kcat = 186.6
    km_ugt1a3_e1_e13g <- fixed(77000)
    label("Michaelis constant, E1 -> E1-3-glucuronide by UGT1A3 (nM)") # Additional file 1 v13 Km; deposit km = 77000
    kcat_sult1e1_e1_e13s <- fixed(150)
    label("Turnover number, E1 -> E1-3-sulfate by SULT1E1 (1/h)") # Additional file 1 v14 kcat; deposit Kcat = 150; kcat: supplement v14 prints 8.46e+6 h-1
    km_sult1e1_e1_e13s <- fixed(6)
    label("Michaelis constant, E1 -> E1-3-sulfate by SULT1E1 (nM)") # Additional file 1 v14 Km; deposit km = 6
    kcat_cyp1a2_e1_oh2e1 <- fixed(552)
    label("Turnover number, E1 -> 2-hydroxy-E1 by CYP1A2 (1/h)") # Additional file 1 v8 kcat; deposit Kcat = 552
    km_cyp1a2_e1_oh2e1 <- fixed(19000)
    label("Michaelis constant, E1 -> 2-hydroxy-E1 by CYP1A2 (nM)") # Additional file 1 v8 Km; deposit km = 19000; Km: supplement v8 prints 1900 nM
    kcat_cyp3a4_e1_oh2e1 <- fixed(42)
    label("Turnover number, E1 -> 2-hydroxy-E1 by CYP3A4 (1/h)") # Additional file 1 v9 kcat; deposit Kcat = 42
    km_cyp3a4_e1_oh2e1 <- fixed(102000)
    label("Michaelis constant, E1 -> 2-hydroxy-E1 by CYP3A4 (nM)") # Additional file 1 v9 Km; deposit km = 102000
    kcat_cyp1a2_e1_oh4e1 <- fixed(120)
    label("Turnover number, E1 -> 4-hydroxy-E1 by CYP1A2 (1/h)") # Additional file 1 v10 kcat; deposit Kcat = 120
    km_cyp1a2_e1_oh4e1 <- fixed(17000)
    label("Michaelis constant, E1 -> 4-hydroxy-E1 by CYP1A2 (nM)") # Additional file 1 v10 Km; deposit km = 17000
    kcat_cyp3a4_e1_oh4e1 <- fixed(30)
    label("Turnover number, E1 -> 4-hydroxy-E1 by CYP3A4 (1/h)") # Additional file 1 v11 kcat; deposit Kcat = 30
    km_cyp3a4_e1_oh4e1 <- fixed(78000)
    label("Michaelis constant, E1 -> 4-hydroxy-E1 by CYP3A4 (nM)") # Additional file 1 v11 Km; deposit km = 78000
    kcat_cyp3a4_e1_oh16ae1 <- fixed(30)
    label("Turnover number, E1 -> 16alpha-hydroxy-E1 by CYP3A4 (1/h)") # Additional file 1 v7 kcat; deposit Kcat = 30
    km_cyp3a4_e1_oh16ae1 <- fixed(64000)
    label("Michaelis constant, E1 -> 16alpha-hydroxy-E1 by CYP3A4 (nM)") # Additional file 1 v7 Km; deposit km = 64000
    kcat_comt_oh2e2_me2e2 <- fixed(408)
    label("Turnover number, 2-hydroxy-E2 -> 2-methoxy-E2 by COMT (1/h)") # Additional file 1 v4 kcat; deposit Kcat = 408
    km_comt_oh2e2_me2e2 <- fixed(108000)
    label("Michaelis constant, 2-hydroxy-E2 -> 2-methoxy-E2 by COMT (nM)") # Additional file 1 v4 Km; deposit km = 108000
    kcat_comt_oh4e2_me4e2 <- fixed(204)
    label("Turnover number, 4-hydroxy-E2 -> 4-methoxy-E2 by COMT (1/h)") # Additional file 1 v6 kcat; deposit Kcat = 204
    km_comt_oh4e2_me4e2 <- fixed(24000)
    label("Michaelis constant, 4-hydroxy-E2 -> 4-methoxy-E2 by COMT (nM)") # Additional file 1 v6 Km; deposit km = 24000
    kcat_comt_oh2e1_me2e1 <- fixed(198)
    label("Turnover number, 2-hydroxy-E1 -> 2-methoxy-E1 by COMT (1/h)") # Additional file 1 v3 kcat; deposit Kcat = 198
    km_comt_oh2e1_me2e1 <- fixed(7400)
    label("Michaelis constant, 2-hydroxy-E1 -> 2-methoxy-E1 by COMT (nM)") # Additional file 1 v3 Km; deposit km = 7400
    kcat_comt_oh4e1_me4e1 <- fixed(402)
    label("Turnover number, 4-hydroxy-E1 -> 4-methoxy-E1 by COMT (1/h)") # Additional file 1 v5 kcat; deposit Kcat = 402
    km_comt_oh4e1_me4e1 <- fixed(53000)
    label("Michaelis constant, 4-hydroxy-E1 -> 4-methoxy-E1 by COMT (nM)") # Additional file 1 v5 Km; deposit km = 53000
    kcat_abcg2_e13s_bile <- fixed(8460000)
    label("Turnover number, E1-3-sulfate transport to bile by ABCG2 (1/h)") # Additional file 1 v106 kcat; deposit Kcat = 8460000; kcat: supplement v106 prints 150 h-1
    km_abcg2_e13s_bile <- fixed(17000)
    label("Michaelis constant, E1-3-sulfate transport to bile by ABCG2 (nM)") # Additional file 1 v106 Km; deposit km = 17000
    kcat_mrp2_e217bg_bile <- fixed(1.9e+07)
    label("Turnover number, E2-17beta-glucuronide transport to bile by MRP2 (1/h)") # Additional file 1 v121 kcat; deposit Kcat = 1.9e+07
    km_mrp2_e217bg_bile <- fixed(98000)
    label("Michaelis constant, E2-17beta-glucuronide transport to bile by MRP2 (nM)") # Additional file 1 v121 Km; deposit km = 98000
    kcat_mdr1_e217bg_bile <- fixed(61981)
    label("Turnover number, E2-17beta-glucuronide transport to bile by MDR1 (1/h)") # Additional file 1 v120 kcat; deposit Kcat = 61981
    km_mdr1_e217bg_bile <- fixed(62000)
    label("Michaelis constant, E2-17beta-glucuronide transport to bile by MDR1 (nM)") # Additional file 1 v120 Km; deposit km = 62000; Km: supplement v120 prints 6200 nM
    kcat_abcg2_e217bg_bile <- fixed(372417)
    label("Turnover number, E2-17beta-glucuronide transport to bile by ABCG2 (1/h)") # Additional file 1 v119 kcat; deposit Kcat = 372417
    km_abcg2_e217bg_bile <- fixed(44200)
    label("Michaelis constant, E2-17beta-glucuronide transport to bile by ABCG2 (nM)") # Additional file 1 v119 Km; deposit km = 44200
    kcat_mrp3_e217bg_blood <- fixed(62976)
    label("Turnover number, E2-17beta-glucuronide transport to liver blood by MRP3 (1/h)") # Additional file 1 v122 kcat; deposit Kcat = 62976
    km_mrp3_e217bg_blood <- fixed(26000)
    label("Michaelis constant, E2-17beta-glucuronide transport to liver blood by MRP3 (nM)") # Additional file 1 v122 Km; deposit km = 26000
    kcat_mrp4_e217bg_blood <- fixed(305035)
    label("Turnover number, E2-17beta-glucuronide transport to liver blood by MRP4 (1/h)") # Additional file 1 v123 kcat; deposit Kcat = 305035
    km_mrp4_e217bg_blood <- fixed(30000)
    label("Michaelis constant, E2-17beta-glucuronide transport to liver blood by MRP4 (nM)") # Additional file 1 v123 Km; deposit km = 30000
  })
  model({
    # ---- Anthropometrics ----
    ht_m <- HT / 100 # source equations use height in metres
    bsa <- 0.007184 * WT^0.425 * HT^0.725 # DuBois BSA (m^2), height in cm
    bmi <- WT / ht_m^2
    # Literal approximations of e exactly as deposited.
    e_a <- 2.718
    e_b <- 2.71

    # ---- Organ volumes, L (Additional file 1 Global Quantities; Bosgra 2012) ----
    v_adipose <- WT * ((1.2 * bmi + 0.23 * AGE - 5.4) / 100) / 0.916
    v_plasma_arterial <- ((2.66 * bsa - 0.46) * 0.33 * 0.55) / 1.025
    v_plasma_venous <- ((2.66 * bsa - 0.46) * 0.67 * 0.55) / 1.025
    v_bone <- e_a^(2.67 * log(ht_m) + 0.0689) / 1.99
    v_brain <- (0.373 * (3.68 - 2.68 * e_b^(-AGE / 0.89)) * (e_b^(-AGE / 629))) / 1.035
    v_gonads <- ((3.3 + 90 * (1 - e_b^(-(AGE / 16.8)^6.7))) / 1000) / 1.03
    v_heart <- e_a^(2.13 * log(ht_m) - 2.502) / 1.03
    v_intestine <- e_a^(2.47 * log(ht_m) - 1.351) / 1.042
    v_kidney <- e_a^(1.93 * log(ht_m) - 2.306) / 1.05
    v_liver <- e_a^(1.98 * log(ht_m) - 0.6786) / 1.03
    v_lung <- e_a^(2.1 * log(ht_m) - 2.092) / 1.05
    v_pancreas <- e_a^(2.43 * log(ht_m) - 3.431) / 1.045
    v_other <- e_a^(1.95 * log(ht_m) - 0.072) / 1.03
    v_skin <- e_b^(1.63 * bsa - 1.93) / 1.5
    v_spleen <- e_a^(2.16 * log(ht_m) - 3.123) / 1.054
    v_stomach <- e_a^(2.45 * log(ht_m) - 3.266) / 1.044
    v_thymus <- ((14 * (7.1 - 6.1 * e_b^(-AGE / 11.9)) * (0.14 + 0.86 * e_b^(-AGE / 10.3))) / 1000) / 1.03
    v_muscle <- (WT * 0.93 - (v_adipose + v_plasma_arterial + v_plasma_venous +
      v_bone + v_brain + v_gonads + v_heart + v_intestine + v_kidney + v_liver +
      v_lung + v_pancreas + v_other + v_skin + v_spleen + v_stomach + v_thymus)) / 1.041

    # ---- Blood flows, L/h: volume * specific flow * 0.55 ----
    # As deposited in MODEL1711210002: unlike the base PBPK variant there is
    # no tissue-density factor, and lung flow is 0.55 * cardiac output.
    q_adipose <- v_adipose * 1.4 * 0.55
    q_bone <- v_bone * 1.8 * 0.55
    q_brain <- v_brain * 32.7 * 0.55
    q_gonads <- v_gonads * 0.8 * 0.55
    q_heart <- v_heart * 70.8 * 0.55
    q_intestine <- v_intestine * 59 * 0.55
    q_kidney <- v_kidney * 215 * 0.55
    q_liver_artery <- v_liver * 16.4 * 0.55
    q_muscle <- v_muscle * 2.4 * 0.55
    q_pancreas <- v_pancreas * 29.5 * 0.55
    q_other <- v_other * 18.8 * 0.55
    q_skin <- v_skin * 7.7 * 0.55
    q_spleen <- v_spleen * 81.5 * 0.55
    q_stomach <- v_stomach * 25.3 * 0.55
    q_thymus <- v_thymus * 18.5 * 0.55
    q_liver_vein <- q_liver_artery + q_intestine + q_pancreas + q_spleen + q_stomach
    q_cardiac <- q_adipose + q_bone + q_brain + q_gonads + q_heart + q_intestine +
      q_kidney + q_liver_artery + q_muscle + q_pancreas + q_other + q_skin +
      q_spleen + q_stomach + q_thymus
    q_lung <- q_cardiac * 0.55

    ka <- exp(lka)

    # ---- Enzyme and transporter amounts (nM) ----
    comt <- abund_comt * f_protein_expression
    cyp1a2 <- abund_cyp1a2 * f_protein_expression
    cyp3a4 <- abund_cyp3a4 * f_protein_expression
    hsd17b1 <- abund_hsd17b1 * f_protein_expression
    hsd17b2 <- abund_hsd17b2 * f_protein_expression
    sult1e1 <- abund_sult1e1 * f_protein_expression
    sult2a1 <- abund_sult2a1 * f_protein_expression
    ugt1a1 <- abund_ugt1a1 * f_protein_expression
    ugt1a3 <- abund_ugt1a3 * f_protein_expression
    ugt2b7 <- abund_ugt2b7 * f_protein_expression
    abcg2 <- abund_abcg2
    mdr1 <- abund_mdr1 * f_protein_expression
    mrp2 <- abund_mrp2 * f_protein_expression
    mrp3 <- abund_mrp3 * f_protein_expression
    mrp4 <- abund_mrp4 * f_protein_expression

    # ---- PBPK fluxes (nM/h), each flow divided by the SOURCE volume as deposited ----
    f_lung_in <- q_lung / v_lung * venous
    f_lung_out <- q_lung / v_lung * lung
    f_adipose_in <- q_adipose / v_adipose * arterial
    f_adipose_out <- q_adipose / v_adipose * adipose
    f_bone_in <- q_bone / v_bone * arterial
    f_bone_out <- q_bone / v_bone * bone
    f_brain_in <- q_brain / v_brain * arterial
    f_brain_out <- q_brain / v_brain * brain
    f_heart_in <- q_heart / v_heart * arterial
    f_heart_out <- q_heart / v_heart * heart
    f_muscle_in <- q_muscle / v_muscle * arterial
    f_muscle_out <- q_muscle / v_muscle * muscle
    f_skin_in <- q_skin / v_skin * arterial
    f_skin_out <- q_skin / v_skin * skin
    f_kidney_in <- q_kidney / v_kidney * arterial
    f_kidney_out <- q_kidney / v_kidney * kidney
    f_thymus_in <- q_thymus / v_thymus * arterial
    f_thymus_out <- q_thymus / v_thymus * thymus
    f_other_in <- q_other / v_other * arterial
    f_other_out <- q_other / v_other * other
    f_gonads_in <- q_gonads / v_gonads * arterial
    f_gonads_out <- q_gonads / v_gonads * vp_gonads
    f_intestine_in <- q_intestine / v_intestine * arterial
    f_intestine_liver <- q_intestine / v_intestine * intestine
    f_spleen_in <- q_spleen / v_spleen * arterial
    f_spleen_liver <- q_spleen / v_spleen * spleen
    f_stomach_in <- q_stomach / v_stomach * arterial
    f_stomach_liver <- q_stomach / v_stomach * stomach
    f_pancreas_in <- q_pancreas / v_pancreas * arterial
    f_pancreas_liver <- q_pancreas / v_pancreas * pancreas
    f_liver_artery <- q_liver_artery / v_liver * arterial
    f_liver_vein <- q_liver_vein / v_liver * vp_liver

    f_perm_liver <- kperm_in_liver * vp_liver - kperm_out_liver * liver
    f_perm_gonads <- kperm_in_gonads * vp_gonads - kperm_out_gonads * gonads
    f_er_liver <- kon_er_liver * target_er_liver * liver - koff_er_liver * complex_er_liver
    f_er_gonads <- kon_er_gonads * target_er_gonads * gonads - koff_er_gonads * complex_er_gonads
    f_nsb_liver <- kon_nsb * liver * nsb_liver - koff_nsb * complex_nsb_liver
    f_shbg_venous <- kon_shbg * venous * shbg_venous - koff_shbg * complex_shbg_venous
    f_albumin_venous <- kon_albumin * venous * albumin_venous - koff_albumin * complex_albumin_venous
    f_shbg_arterial <- kon_shbg * arterial * shbg_arterial - koff_shbg * complex_shbg_arterial
    f_albumin_arterial <- kon_albumin * arterial * albumin_arterial - koff_albumin * complex_albumin_arterial
    f_cleh <- kcleh * kidney
    f_iv <- kiv * depot_iv
    f_oral <- ka * depot

    # ---- Liver metabolism and transport (nM/h); vNN = Additional file 1 LiverODE reaction ----
    r_ugt1a1_e2_e23g <- liver * ugt1a1 * kcat_ugt1a1_e2_e23g / (liver + km_ugt1a1_e2_e23g) # v27
    r_ugt1a3_e2_e23g <- liver * ugt1a3 * kcat_ugt1a3_e2_e23g / (liver + km_ugt1a3_e2_e23g) # v28
    r_ugt1a1_e2_e217bg <- liver * ugt1a1 * kcat_ugt1a1_e2_e217bg / (liver + km_ugt1a1_e2_e217bg) # v24
    r_ugt1a3_e2_e217bg <- liver * ugt1a3 * kcat_ugt1a3_e2_e217bg / (liver + km_ugt1a3_e2_e217bg) # v25
    r_ugt2b7_e2_e217bg <- liver * ugt2b7 * kcat_ugt2b7_e2_e217bg / (liver + km_ugt2b7_e2_e217bg) # v26
    r_sult2a1_e2_e23s <- liver * sult2a1 * kcat_sult2a1_e2_e23s / (liver + km_sult2a1_e2_e23s) # v29
    r_cyp1a2_e2_oh2e2 <- liver * cyp1a2 * kcat_cyp1a2_e2_oh2e2 / (liver + km_cyp1a2_e2_oh2e2) # v18
    r_cyp3a4_e2_oh2e2 <- liver * cyp3a4 * kcat_cyp3a4_e2_oh2e2 / (liver + km_cyp3a4_e2_oh2e2) # v19
    r_cyp1a2_e2_oh4e2 <- liver * cyp1a2 * kcat_cyp1a2_e2_oh4e2 / (liver + km_cyp1a2_e2_oh4e2) # v20
    r_cyp3a4_e2_oh4e2 <- liver * cyp3a4 * kcat_cyp3a4_e2_oh4e2 / (liver + km_cyp3a4_e2_oh4e2) # v21
    r_cyp1a2_e2_oh16ae2 <- liver * cyp1a2 * kcat_cyp1a2_e2_oh16ae2 / (liver + km_cyp1a2_e2_oh16ae2) # v16
    r_cyp3a4_e2_oh16ae2 <- liver * cyp3a4 * kcat_cyp3a4_e2_oh16ae2 / (liver + km_cyp3a4_e2_oh16ae2) # v17
    r_hsd17b1_e2_e1 <- liver * hsd17b1 * kcat_hsd17b1_e2_e1 / (liver + km_hsd17b1_e2_e1) # v22
    r_hsd17b2_e2_e1 <- liver * hsd17b2 * kcat_hsd17b2_e2_e1 / (liver + km_hsd17b2_e2_e1) # v23
    r_hsd17b1_e1_e2 <- e1_liver * hsd17b1 * kcat_hsd17b1_e1_e2 / (e1_liver + km_hsd17b1_e1_e2) # v15
    r_ugt1a1_e1_e13g <- e1_liver * ugt1a1 * kcat_ugt1a1_e1_e13g / (e1_liver + km_ugt1a1_e1_e13g) # v12
    r_ugt1a3_e1_e13g <- e1_liver * ugt1a3 * kcat_ugt1a3_e1_e13g / (e1_liver + km_ugt1a3_e1_e13g) # v13
    r_sult1e1_e1_e13s <- e1_liver * sult1e1 * kcat_sult1e1_e1_e13s / (e1_liver + km_sult1e1_e1_e13s) # v14
    r_cyp1a2_e1_oh2e1 <- e1_liver * cyp1a2 * kcat_cyp1a2_e1_oh2e1 / (e1_liver + km_cyp1a2_e1_oh2e1) # v8
    r_cyp3a4_e1_oh2e1 <- e1_liver * cyp3a4 * kcat_cyp3a4_e1_oh2e1 / (e1_liver + km_cyp3a4_e1_oh2e1) # v9
    r_cyp1a2_e1_oh4e1 <- e1_liver * cyp1a2 * kcat_cyp1a2_e1_oh4e1 / (e1_liver + km_cyp1a2_e1_oh4e1) # v10
    r_cyp3a4_e1_oh4e1 <- e1_liver * cyp3a4 * kcat_cyp3a4_e1_oh4e1 / (e1_liver + km_cyp3a4_e1_oh4e1) # v11
    r_cyp3a4_e1_oh16ae1 <- e1_liver * cyp3a4 * kcat_cyp3a4_e1_oh16ae1 / (e1_liver + km_cyp3a4_e1_oh16ae1) # v7
    r_comt_oh2e2_me2e2 <- oh2e2_liver * comt * kcat_comt_oh2e2_me2e2 / (oh2e2_liver + km_comt_oh2e2_me2e2) # v4
    r_comt_oh4e2_me4e2 <- oh4e2_liver * comt * kcat_comt_oh4e2_me4e2 / (oh4e2_liver + km_comt_oh4e2_me4e2) # v6
    r_comt_oh2e1_me2e1 <- oh2e1_liver * comt * kcat_comt_oh2e1_me2e1 / (oh2e1_liver + km_comt_oh2e1_me2e1) # v3
    r_comt_oh4e1_me4e1 <- oh4e1_liver * comt * kcat_comt_oh4e1_me4e1 / (oh4e1_liver + km_comt_oh4e1_me4e1) # v5
    r_abcg2_e13s_bile <- e13s_liver * abcg2 * kcat_abcg2_e13s_bile / (e13s_liver + km_abcg2_e13s_bile) # v106
    r_mrp2_e217bg_bile <- e217bg_liver * mrp2 * kcat_mrp2_e217bg_bile / (e217bg_liver + km_mrp2_e217bg_bile) # v121
    r_mdr1_e217bg_bile <- e217bg_liver * mdr1 * kcat_mdr1_e217bg_bile / (e217bg_liver + km_mdr1_e217bg_bile) # v120
    r_abcg2_e217bg_bile <- e217bg_liver * abcg2 * kcat_abcg2_e217bg_bile / (e217bg_liver + km_abcg2_e217bg_bile) # v119
    r_mrp3_e217bg_blood <- e217bg_liver * mrp3 * kcat_mrp3_e217bg_blood / (e217bg_liver + km_mrp3_e217bg_blood) # v122
    r_mrp4_e217bg_blood <- e217bg_liver * mrp4 * kcat_mrp4_e217bg_blood / (e217bg_liver + km_mrp4_e217bg_blood) # v123

    # ---- Initial conditions: deposited endogenous steady state ----
    # (BioModels MODEL1711210002 initial concentrations, reference woman
    # 60 kg / 1.63 m / 35 y.)
    venous(0) <- 0.0209437464725829
    complex_shbg_venous(0) <- 0.550809233315276
    shbg_venous(0) <- bl_shbg - 0.550809233315276
    complex_albumin_venous(0) <- 0.615991785586703
    albumin_venous(0) <- bl_albumin - 0.615991785586703
    arterial(0) <- 0.020321191153055
    complex_shbg_arterial(0) <- 0.534655210167046
    shbg_arterial(0) <- bl_shbg - 0.534655210167046
    complex_albumin_arterial(0) <- 0.597681378290001
    albumin_arterial(0) <- bl_albumin - 0.597681378290001
    lung(0) <- 0.0209437464725832
    adipose(0) <- 0.020321191153055
    bone(0) <- 0.020321191153055
    brain(0) <- 0.020321191153055
    heart(0) <- 0.020321191153055
    muscle(0) <- 0.020321191153055
    skin(0) <- 0.020321191153055
    intestine(0) <- 0.020321191153055
    spleen(0) <- 0.020321191153055
    stomach(0) <- 0.020321191153055
    pancreas(0) <- 0.020321191153055
    thymus(0) <- 0.020321191153055
    other(0) <- 0.020321191153055
    kidney(0) <- 0.0194968020596256
    vp_liver(0) <- 0.0126846706056808
    liver(0) <- 4.42687933551509e-06
    complex_er_liver(0) <- 0.00081453137446636
    target_er_liver(0) <- bl_er_liver - 0.00081453137446636
    complex_nsb_liver(0) <- 2.51633094477505e-06
    nsb_liver(0) <- bl_nsb_liver - 2.51633094477505e-06
    e1_liver(0) <- 1.39373733418986e-08
    oh2e2_liver(0) <- 0.000918301772841599
    oh4e2_liver(0) <- 2.48361819955745e-05
    oh2e1_liver(0) <- 3.56059427187798e-07
    oh4e1_liver(0) <- 3.14143241763599e-07
    e13s_liver(0) <- 2.47334883947458e-08
    e217bg_liver(0) <- 2.49989643894842e-08
    vp_gonads(0) <- 4.56577573721945
    gonads(0) <- 45.677757372197
    complex_er_gonads(0) <- 139.237933
    target_er_gonads(0) <- bl_er_gonads - 139.237933

    # ---- Balance equations ----
    d/dt(depot) <- -f_oral - kfeces * depot
    d/dt(depot_iv) <- -f_iv
    d/dt(venous) <- f_adipose_out + f_bone_out + f_brain_out + f_heart_out +
      f_muscle_out + f_skin_out + f_liver_vein + f_kidney_out + f_thymus_out +
      f_gonads_out + f_other_out + f_iv - f_lung_in - f_shbg_venous - f_albumin_venous
    d/dt(shbg_venous) <- -f_shbg_venous
    d/dt(complex_shbg_venous) <- f_shbg_venous
    d/dt(albumin_venous) <- -f_albumin_venous
    d/dt(complex_albumin_venous) <- f_albumin_venous
    d/dt(arterial) <- f_lung_out - f_adipose_in - f_bone_in - f_brain_in -
      f_heart_in - f_muscle_in - f_skin_in - f_intestine_in - f_liver_artery -
      f_spleen_in - f_kidney_in - f_thymus_in - f_gonads_in - f_stomach_in -
      f_pancreas_in - f_other_in - f_shbg_arterial - f_albumin_arterial
    d/dt(shbg_arterial) <- -f_shbg_arterial
    d/dt(complex_shbg_arterial) <- f_shbg_arterial
    d/dt(albumin_arterial) <- -f_albumin_arterial
    d/dt(complex_albumin_arterial) <- f_albumin_arterial
    d/dt(lung) <- f_lung_in - f_lung_out
    d/dt(adipose) <- f_adipose_in - f_adipose_out
    d/dt(bone) <- f_bone_in - f_bone_out
    d/dt(brain) <- f_brain_in - f_brain_out
    d/dt(heart) <- f_heart_in - f_heart_out
    d/dt(muscle) <- f_muscle_in - f_muscle_out
    d/dt(skin) <- f_skin_in - f_skin_out
    d/dt(intestine) <- f_intestine_in + f_oral - f_intestine_liver
    d/dt(spleen) <- f_spleen_in - f_spleen_liver
    d/dt(stomach) <- f_stomach_in - f_stomach_liver
    d/dt(pancreas) <- f_pancreas_in - f_pancreas_liver
    d/dt(thymus) <- f_thymus_in - f_thymus_out
    d/dt(kidney) <- f_kidney_in - f_kidney_out - f_cleh
    d/dt(other) <- f_other_in - f_other_out
    d/dt(vp_liver) <- f_liver_artery + f_intestine_liver + f_spleen_liver +
      f_stomach_liver + f_pancreas_liver - f_liver_vein - f_perm_liver
    # Hepatocyte E2: permeability uptake, ER and non-specific binding, the
    # E1 <-> E2 interconversion and every E2-consuming conversion.
    d/dt(liver) <- f_perm_liver - f_er_liver - f_nsb_liver +
      r_hsd17b1_e1_e2 - r_hsd17b1_e2_e1 - r_hsd17b2_e2_e1 -
      r_ugt1a1_e2_e23g - r_ugt1a3_e2_e23g -
      r_ugt1a1_e2_e217bg - r_ugt1a3_e2_e217bg - r_ugt2b7_e2_e217bg -
      r_sult2a1_e2_e23s -
      r_cyp1a2_e2_oh2e2 - r_cyp3a4_e2_oh2e2 -
      r_cyp1a2_e2_oh4e2 - r_cyp3a4_e2_oh4e2 -
      r_cyp1a2_e2_oh16ae2 - r_cyp3a4_e2_oh16ae2
    d/dt(target_er_liver) <- -f_er_liver
    d/dt(complex_er_liver) <- f_er_liver
    d/dt(nsb_liver) <- -f_nsb_liver
    d/dt(complex_nsb_liver) <- f_nsb_liver
    d/dt(e1_liver) <- r_hsd17b1_e2_e1 + r_hsd17b2_e2_e1 - r_hsd17b1_e1_e2 -
      r_ugt1a1_e1_e13g - r_ugt1a3_e1_e13g - r_sult1e1_e1_e13s -
      r_cyp1a2_e1_oh2e1 - r_cyp3a4_e1_oh2e1 -
      r_cyp1a2_e1_oh4e1 - r_cyp3a4_e1_oh4e1 - r_cyp3a4_e1_oh16ae1
    d/dt(oh2e2_liver) <- r_cyp1a2_e2_oh2e2 + r_cyp3a4_e2_oh2e2 - r_comt_oh2e2_me2e2
    d/dt(oh4e2_liver) <- r_cyp1a2_e2_oh4e2 + r_cyp3a4_e2_oh4e2 - r_comt_oh4e2_me4e2
    d/dt(oh2e1_liver) <- r_cyp1a2_e1_oh2e1 + r_cyp3a4_e1_oh2e1 - r_comt_oh2e1_me2e1
    d/dt(oh4e1_liver) <- r_cyp1a2_e1_oh4e1 + r_cyp3a4_e1_oh4e1 - r_comt_oh4e1_me4e1
    d/dt(e13s_liver) <- r_sult1e1_e1_e13s - r_abcg2_e13s_bile
    d/dt(e217bg_liver) <- r_ugt1a1_e2_e217bg + r_ugt1a3_e2_e217bg + r_ugt2b7_e2_e217bg -
      r_mrp2_e217bg_bile - r_mdr1_e217bg_bile - r_abcg2_e217bg_bile -
      r_mrp3_e217bg_blood - r_mrp4_e217bg_blood
    d/dt(vp_gonads) <- f_gonads_in - f_gonads_out - f_perm_gonads
    d/dt(gonads) <- f_perm_gonads - f_er_gonads + ksyn_gonads
    d/dt(target_er_gonads) <- -f_er_gonads
    d/dt(complex_er_gonads) <- f_er_gonads

    # Free (unbound) venous E2, nM.
    Cc <- venous
  })
}
