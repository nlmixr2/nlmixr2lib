Benson_2014_PF04457845_qsp <- function() {
  description <- paste(
    "QSP. Whole-body four-physiological-compartment (brain, plasma,",
    "rest-of-body ROB, blood-brain-barrier microvascular endothelial",
    "cells MEC) systems-pharmacology model for the irreversible FAAH-1",
    "inhibitor PF-04457845 and five fatty acid ethanolamide substrates",
    "(anandamide AEA, N-oleoyl-ethanolamide OEA, N-palmitoyl-ethanolamide",
    "PEA, N-linoleoyl-ethanolamide LEA, N-stearoyl-ethanolamide SEA).",
    "39 ODE states coupling: PF-04457845 pharmacokinetics",
    "(first-order absorption, saturable Michaelis-Menten plus linear",
    "elimination, plasma-ROB tissue distribution); irreversible",
    "PF-04457845 inhibition of FAAH with FAAH protein turnover;",
    "NAPE-phospholipase-D-mediated ethanolamide synthesis from the",
    "N-acyl phosphatidylethanolamine precursors with competitive",
    "substrate inhibition among the five ethanolamides; FAAH-independent",
    "clearance via N-acylethanolamine-hydrolyzing acid amidase (NAAA);",
    "bidirectional AEA-transporter and passive diffusion transport of",
    "ethanolamides between ROB, MEC and plasma; blood-brain-barrier",
    "flux of ethanolamides between brain and MEC; and Emax cannabinoid",
    "receptor CB1 occupancy in brain by AEA. All 110+ structural,",
    "physiologic, enzyme-kinetic, receptor-binding, and PK parameters",
    "are fixed from the paper's supplement Table S1 and Model S1 SBML",
    "export: Peters and Hultin 2008 PBPK physiology for the average",
    "70-kg-adult tissue volumes; literature in vitro enzyme-kinetic",
    "measurements for FAAH kcat/Km, NAPE-PLD forward/reverse rates and",
    "substrate Ki, and NAAA affinities; McPartland 2007 CB1 in vitro",
    "K_d; and Pfizer Phase I data (Li 2012) for FAAH turnover, drug PK,",
    "and the fitted NAT relative activities. No IIV, no residual",
    "variability -- the model is deterministic and describes mean",
    "biomarker profiles. 2-arachidonoyl-glycerol (2-AG) is not",
    "included (fixed at zero) per the paper's stated scope. Model was",
    "originally implemented in DBSolve Optimum and Matlab/SimBiology.",
    sep = " "
  )
  reference <- paste(
    "Benson N, Metelkin E, Demin O, Li GL, Nichols D, van der Graaf PH",
    "(2014). A systems pharmacology perspective on the clinical",
    "development of fatty acid amide hydrolase inhibitors for pain.",
    "CPT: Pharmacometrics & Systems Pharmacology 3, e91.",
    "doi:10.1038/psp.2013.72.",
    "Clinical PK/PD data source: Li GL et al. (2012).",
    "Assessment of the pharmacology and tolerability of PF-04457845,",
    "an irreversible inhibitor of fatty acid amide hydrolase-1, in",
    "healthy subjects. Br J Clin Pharmacol 73, 706-716.",
    "PBPK physiology source: Peters SA, Hultin L (2008).",
    "Early identification of drug-induced impairment of gastric",
    "emptying through PBPK simulation. J Pharmacokinet Pharmacodyn 35,",
    "1-30.",
    sep = " "
  )
  vignette <- "Benson_2014_PF04457845_qsp"

  paper_specific_compartments <- c(
    # Brain states (12)
    "aea_brain", "oea_brain", "pea_brain", "lea_brain", "sea_brain",
    "nape_brain", "nope_brain", "nppe_brain", "nlpe_brain", "nspe_brain",
    "faah_brain", "faahinh_brain",
    # Rest-of-body states (12)
    "aea_rob", "oea_rob", "pea_rob", "lea_rob", "sea_rob",
    "nape_rob", "nope_rob", "nppe_rob", "nlpe_rob", "nspe_rob",
    "faah_rob", "faahinh_rob",
    # BBB microvascular endothelial cell states (7)
    "aea_mec", "oea_mec", "pea_mec", "lea_mec", "sea_mec",
    "faah_mec", "faahinh_mec",
    # Plasma ethanolamide states (5)
    "aea_plasma", "oea_plasma", "pea_plasma", "lea_plasma", "sea_plasma",
    # PF-04457845 drug states (3)
    "pf_gut", "pf_p_amt", "pf_r_amt"
  )

  units <- list(
    time          = "h",
    dosing        = paste(
      "PF-04457845 dose into the pf_gut compartment must be in ng",
      "(convert from mg via amt_ng = amt_mg * 1e6).",
      "Bioavailability (Benson 2014 supplement Table S1 saturable Emax:",
      "F_PFM = Emax_PFM * dose_mg / (ED50 + dose_mg)) is applied via",
      "f(pf_gut) in model() using the dose_mg covariate.",
      "The paper's Figure 4 simulated 1 mg and 10 mg single oral doses;",
      "Figure 5 simulated 0.1-40 mg.",
      sep = " "
    ),
    concentration = paste(
      "PF-04457845 plasma concentration Cc is in ng/mL",
      "(via 1e-3 * pf_p_amt / vss_pf where pf_p_amt is amount in ng and",
      "vss_pf is 58.328 L; see Errata for the supplement Table S1 unit",
      "typo 'mL' that the model math requires to be L). Ethanolamide",
      "plasma concentrations Cc_AEA/Cc_OEA/Cc_PEA/Cc_LEA/Cc_SEA are in",
      "nmol/L (nM). Ethanolamide and NAPE intra-tissue concentrations",
      "are in nmol/L. FAAH and FAAHinh concentrations are in nmol/L.",
      "FAAHact is the plasma-FAAH activity as a fraction of baseline",
      "(0-1). CB1occ is the brain CB1 receptor occupancy fraction (0-1).",
      sep = " "
    )
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. analyte/specimen proposed by a local model from the
  # model description; units derived from the units block. verified = FALSE
  # means NOT checked against the source paper.
  compartmentData <- list(
    aea_brain     = list(analyte = "anandamide AEA", units = NA_character_, specimen = "tissue", verified = FALSE),
    oea_brain     = list(analyte = "N-oleoyl-ethanolamide OEA", units = NA_character_, specimen = "tissue", verified = FALSE),
    pea_brain     = list(analyte = "N-palmitoyl-ethanolamide PEA", units = NA_character_, specimen = "tissue", verified = FALSE),
    lea_brain     = list(analyte = "N-linoleoyl-ethanolamide LEA", units = NA_character_, specimen = "tissue", verified = FALSE),
    sea_brain     = list(analyte = "N-stearoyl-ethanolamide SEA", units = NA_character_, specimen = "tissue", verified = FALSE),
    nape_brain    = list(analyte = "NAPE", units = NA_character_, specimen = "tissue", verified = FALSE),
    nope_brain    = list(analyte = "N-oleoyl-phosphatidylethanolamine NOPE", units = NA_character_, specimen = "tissue", verified = FALSE),
    nppe_brain    = list(analyte = "N-palmitoyl-phosphatidylethanolamine NPPE", units = NA_character_, specimen = "tissue", verified = FALSE),
    nlpe_brain    = list(analyte = "N-linoleoyl-phosphatidylethanolamine NLPE", units = NA_character_, specimen = "tissue", verified = FALSE),
    nspe_brain    = list(analyte = "N-stearoyl-phosphatidylethanolamine NSPE", units = NA_character_, specimen = "tissue", verified = FALSE),
    faah_brain    = list(analyte = "FAAH protein", units = NA_character_, specimen = "tissue", verified = FALSE),
    faahinh_brain = list(analyte = "Inhibited FAAH protein", units = NA_character_, specimen = "tissue", verified = FALSE),
    aea_rob       = list(analyte = "anandamide AEA", units = NA_character_, specimen = "tissue", verified = FALSE),
    oea_rob       = list(analyte = "N-oleoyl-ethanolamide OEA", units = NA_character_, specimen = "tissue", verified = FALSE),
    pea_rob       = list(analyte = "N-palmitoyl-ethanolamide PEA", units = NA_character_, specimen = "tissue", verified = FALSE),
    lea_rob       = list(analyte = "N-linoleoyl-ethanolamide LEA", units = NA_character_, specimen = "tissue", verified = FALSE),
    sea_rob       = list(analyte = "N-stearoyl-ethanolamide SEA", units = NA_character_, specimen = "tissue", verified = FALSE),
    nape_rob      = list(analyte = "NAPE", units = NA_character_, specimen = "tissue", verified = FALSE),
    nope_rob      = list(analyte = "N-oleoyl-phosphatidylethanolamine NOPE", units = NA_character_, specimen = "tissue", verified = FALSE),
    nppe_rob      = list(analyte = "N-palmitoyl-phosphatidylethanolamine NPPE", units = NA_character_, specimen = "tissue", verified = FALSE),
    nlpe_rob      = list(analyte = "N-linoleoyl-phosphatidylethanolamine NLPE", units = NA_character_, specimen = "tissue", verified = FALSE),
    nspe_rob      = list(analyte = "N-stearoyl-phosphatidylethanolamine NSPE", units = NA_character_, specimen = "tissue", verified = FALSE),
    faah_rob      = list(analyte = "FAAH protein", units = NA_character_, specimen = "tissue", verified = FALSE),
    faahinh_rob   = list(analyte = "Inhibited FAAH protein", units = NA_character_, specimen = "tissue", verified = FALSE),
    aea_mec       = list(analyte = "anandamide AEA", units = NA_character_, specimen = "blood cell", verified = FALSE),
    oea_mec       = list(analyte = "N-oleoyl-ethanolamide OEA", units = NA_character_, specimen = "blood cell", verified = FALSE),
    pea_mec       = list(analyte = "N-palmitoyl-ethanolamide PEA", units = NA_character_, specimen = "blood cell", verified = FALSE),
    lea_mec       = list(analyte = "N-linoleoyl-ethanolamide LEA", units = NA_character_, specimen = "blood cell", verified = FALSE),
    sea_mec       = list(analyte = "N-stearoyl-ethanolamide SEA", units = NA_character_, specimen = "blood cell", verified = FALSE),
    faah_mec      = list(analyte = "FAAH protein", units = NA_character_, specimen = "blood cell", verified = FALSE),
    faahinh_mec   = list(analyte = "Inhibited FAAH protein", units = NA_character_, specimen = "blood cell", verified = FALSE),
    aea_plasma    = list(analyte = "anandamide AEA", units = NA_character_, specimen = "plasma", verified = FALSE),
    oea_plasma    = list(analyte = "N-oleoyl-ethanolamide OEA", units = NA_character_, specimen = "plasma", verified = FALSE),
    pea_plasma    = list(analyte = "N-palmitoyl-ethanolamide PEA", units = NA_character_, specimen = "plasma", verified = FALSE),
    lea_plasma    = list(analyte = "N-linoleoyl-ethanolamide LEA", units = NA_character_, specimen = "plasma", verified = FALSE),
    sea_plasma    = list(analyte = "N-stearoyl-ethanolamide SEA", units = NA_character_, specimen = "plasma", verified = FALSE),
    pf_gut        = list(analyte = "PF-04457845", units = NA_character_, specimen = "tissue", verified = FALSE),
    pf_p_amt      = list(analyte = "PF-04457845", units = NA_character_, specimen = "plasma", verified = FALSE),
    pf_r_amt      = list(analyte = "PF-04457845", units = NA_character_, specimen = "tissue", verified = FALSE)
  )

  covariateData <- list(
    DOSE = list(
      description        = paste(
        "PF-04457845 single-dose amount used to compute the",
        "saturable oral bioavailability F_PFM =",
        "Emax_PFM * DOSE / (ED50 + DOSE) (Benson 2014",
        "supplement Table S1). Must be recorded per subject",
        "and per dose event because F_PFM is dose-dependent."
      ),
      units              = "mg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Use case (a) 'per-subject assigned dose' from the DOSE canonical.",
        "The paper fitted the model to single oral doses of 1 mg",
        "and 10 mg (Figure 4). Figure 5 additionally simulated",
        "the 0.1-40 mg range. For a 10 mg dose F_PFM = 0.773 *",
        "10 / (0.53 + 10) = 0.734; for a 1 mg dose F_PFM = 0.773 *",
        "1 / (0.53 + 1) = 0.505."
      ),
      source_name        = "Dose"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = NA_integer_,
    n_studies      = 1L,
    age_range      = "healthy adult",
    weight_range   = "average 70 kg (Peters & Hultin 2008 reference PBPK physiology)",
    sex_female_pct = NA_real_,
    race_ethnicity = NA_character_,
    disease_state  = paste(
      "Healthy adult subjects in the PF-04457845 Phase I clinical",
      "trial reported by Li et al. 2012 (Br J Clin Pharmacol",
      "73, 706-716). PF-04457845 was originally being developed",
      "for the treatment of pain (subsequently evaluated in",
      "osteoarthritis Phase II by Huggins et al. 2012 Pain 153,",
      "1837-1846 with no analgesic effect observed)."
    ),
    dose_range     = paste(
      "Single oral doses of 1 mg and 10 mg were fitted (Benson",
      "2014 Figure 4). The paper additionally simulated CB1",
      "occupancy at doses of 0.1, 1, 5, 10, 20, and 40 mg",
      "(Figure 5)."
    ),
    regions        = NA_character_,
    notes          = paste(
      "The four physiological compartment volumes (BRAIN = 1.45 L,",
      "PLASMA = 2.649 L, ROB = 65.3 L, MEC = 1.5e-5 L) and the twelve",
      "individual tissue sub-compartment volumes (LIVER, Gut, Spleen,",
      "Kidney, Lungs, Heart, Muscles, Pancreas, Testis, Thymus,",
      "Leucocytes, Brain) are taken from Peters & Hultin 2008 for an",
      "average 70-kg adult. All 110+ parameter values are FIXED; no",
      "IIV or residual variability are estimated by the paper. The",
      "2-arachidonoyl-glycerol brain concentration ag2_brain (which",
      "would competitively occupy CB1 alongside AEA) and its binding",
      "constant kd_ag2 are fixed at zero per the paper's stated model",
      "scope (Section on Discussion, 'A limitation of our model is",
      "that it does not include the available data on the CB agonist",
      "2-arachidonoyl glycerol'). See Errata for six numeric",
      "supplement-vs-SBML discrepancies where the executable SBML",
      "Model S1 values were used in preference to the supplement",
      "Table S1 values."
    )
  )

  ini({
    # =====================================================================
    # ALL 110+ parameters are FIXED. Values from Benson 2014 supplement
    # Model S1 SBML export (where the SBML and Table S1 disagree, the
    # executable SBML value is used and the discrepancy is documented in
    # the vignette Errata section). Sources follow Table S1 "Basis of the
    # choice" / "Source of information" columns.
    # =====================================================================

    # ---------- Physiological compartment volumes (L) --------------------
    brain_vol      <- fixed(1.45)        ; label("Brain volume (L; Peters 2008 PBPK)")                    # Table S1 BRAIN = 1.45 L; Peters 2008 J Pharmacokinet Pharmacodyn 35:1-30
    plasma_vol     <- fixed(2.649)       ; label("Plasma volume (L; Peters 2008 PBPK)")                   # Table S1 PLASMA = 2.649 L
    rob_vol        <- fixed(65.3)        ; label("Rest-of-body volume (L; Peters 2008 PBPK)")             # Table S1 ROB = 65.3 L (= 70 kg total minus brain minus plasma minus MEC)
    mec_vol        <- fixed(1.5e-5)      ; label("BBB microvascular endothelial cell volume (L)")         # Table S1 MEC = 1.5e-5 L (calculated)

    # ---------- Individual tissue sub-volumes (L) that feed c_NAT_ROB,
    # ---------- c_FAAH_ROB and c_NAAA_ROB weighted sums (Peters 2008) ----
    liver_vol      <- fixed(1.69)        ; label("Liver volume (L; SBML Model S1)")                       # SBML Model S1 LIVER = 1.69 L (Errata: Table S1 lists 65.3 L which duplicates ROB and is a supplement typo -- SBML value used)
    gut_vol        <- fixed(1.65)        ; label("Gut volume (L; Peters 2008 PBPK)")                      # Table S1 Gut = 1.65 L
    spleen_vol     <- fixed(0.192)       ; label("Spleen volume (L; Peters 2008 PBPK)")                   # Table S1 Spleen = 0.192 L
    kidney_vol     <- fixed(0.28)        ; label("Kidney volume (L; Peters 2008 PBPK)")                   # Table S1 Kidney = 0.28 L
    lungs_vol      <- fixed(1.172)       ; label("Lungs volume (L; Peters 2008 PBPK)")                    # Table S1 Lungs = 1.172 L
    heart_vol      <- fixed(0.31)        ; label("Heart volume (L; Peters 2008 PBPK)")                    # Table S1 Heart = 0.31 L
    muscles_vol    <- fixed(35.0)        ; label("Muscles volume (L; Peters 2008 PBPK)")                  # Table S1 Muscles = 35 L
    pancreas_vol   <- fixed(0.077)       ; label("Pancreas volume (L; Peters 2008 PBPK)")                 # Table S1 Pancreas = 0.077 L
    testis_vol     <- fixed(0.036)       ; label("Testis volume (L; Peters 2008 PBPK)")                   # Table S1 Testis = 0.036 L
    thymus_vol     <- fixed(0.029)       ; label("Thymus volume (L; Peters 2008 PBPK)")                   # Table S1 Thymus = 0.029 L
    leucocytes_vol <- fixed(0.025)       ; label("Leucocytes volume (L; calculated blood/leucocyte ratio)") # Table S1 Leucocytes = 0.025 L (blood/leucocyte volume ratio 200; CMAJ 1967;97:793-6)

    # ---------- NAPE synthesis by N-acyltransferase (NAT), all
    # ---------- ethanolamides ------------------------------------------
    vmax_nat       <- fixed(300)         ; label("NAT maximum rate Vmax_NAT (nmol/L/h; fitted to Phase I data)")   # Table S1 Vmax_NAT = 300 nM/h; fitted Pfizer CT I
    a_nat_a        <- fixed(1)           ; label("NAT relative activity toward AEA precursor a_NAT_A (unitless)")  # Table S1 a_NAT_A = 1 (fixed reference)
    a_nat_o        <- fixed(13)          ; label("NAT relative activity toward OEA precursor a_NAT_O (unitless)")  # Table S1 a_NAT_O = 13; fitted Pfizer CT I
    a_nat_p        <- fixed(0.42)        ; label("NAT relative activity toward PEA precursor a_NAT_P (unitless)")  # Table S1 a_NAT_P = 0.42; fitted Pfizer CT I
    a_nat_l        <- fixed(8.6)         ; label("NAT relative activity toward LEA precursor a_NAT_L (unitless)")  # Table S1 a_NAT_L = 8.6; fitted Pfizer CT I
    a_nat_s        <- fixed(1)           ; label("NAT relative activity toward SEA precursor a_NAT_S (unitless)")  # Table S1 a_NAT_S = 1 (free / assumed reference)

    # ---------- Fatty acid precursor concentrations (rat testis 1-pos.) --
    p_a            <- fixed(0.051)       ; label("AEA precursor concentration p_A (unitless fraction)")     # Table S1 p_A = 0.051; BBRC 1996 218:113-117
    p_o            <- fixed(0.098)       ; label("OEA precursor concentration p_O (unitless fraction)")     # Table S1 p_O = 0.098; BBRC 1996 218:113-117
    p_p            <- fixed(0.615)       ; label("PEA precursor concentration p_P (unitless fraction; SBML)") # SBML Model S1 p_P = 0.615 (Table S1 rounded to 0.61)
    p_l            <- fixed(0.016)       ; label("LEA precursor concentration p_L (unitless fraction; SBML)") # SBML Model S1 p_L = 0.016 (Table S1 rounded to 0.015)
    p_s            <- fixed(0.191)       ; label("SEA precursor concentration p_S (unitless fraction)")     # Table S1 p_S = 0.191; BBRC 1996 218:113-117

    # ---------- NAT tissue distribution factors b_NAT_<tissue> ----------
    b_nat_brain      <- fixed(1.667)     ; label("NAT tissue factor brain b_NAT_Brain (unitless)")           # Table S1 b_NAT_Brain = 1.667; calculated per J Neurosci 17(4):1226-1242
    b_nat_pancreas   <- fixed(0.333)     ; label("NAT tissue factor pancreas b_NAT_Pancreas (unitless)")     # Table S1 b_NAT_Pancreas = 0.333
    b_nat_kidney     <- fixed(0.667)     ; label("NAT tissue factor kidney b_NAT_Kidney (unitless)")         # Table S1 b_NAT_Kidney = 0.667
    b_nat_heart      <- fixed(1.0)       ; label("NAT tissue factor heart b_NAT_Heart (unitless)")           # Table S1 b_NAT_Heart = 1 (reference)
    b_nat_lungs      <- fixed(0.033)     ; label("NAT tissue factor lungs b_NAT_Lungs (unitless)")           # Table S1 b_NAT_Lungs = 0.033
    b_nat_muscles    <- fixed(0.333)     ; label("NAT tissue factor muscles b_NAT_Muscles (unitless)")       # Table S1 b_NAT_Muscles = 0.333
    b_nat_testis     <- fixed(0.667)     ; label("NAT tissue factor testis b_NAT_Testis (unitless)")         # Table S1 b_NAT_Testis = 0.667
    b_nat_leucocytes <- fixed(0)         ; label("NAT tissue factor leucocytes b_NAT_Leucocytes (unitless)") # SBML Model S1 b_NAT_Leucocytes = 0 (not tabulated in Table S1)

    # ---------- NAPE-PLD forward rate constants and Michaelis constants -
    k_na_pe        <- fixed(202)         ; label("PLD forward rate k_NA_PE for NAPE (1/h)")                  # Table S1 k_NA_PE = 202/h; JBC 2004 279:5298-5305 (in vitro mouse NAPE-PLD)
    k_no_pe        <- fixed(230)         ; label("PLD forward rate k_NO_PE for NOPE (1/h; SBML value)")      # SBML Model S1 k_NO_PE = 230/h (Table S1 rounded to 230.4)
    k_np_pe        <- fixed(270)         ; label("PLD forward rate k_NP_PE for NPPE (1/h)")                  # Table S1 k_NP_PE = 270/h; JBC 2004 279:5298-5305
    k_nl_pe        <- fixed(100)         ; label("PLD forward rate k_NL_PE for NLPE (1/h; free/assumed)")    # Table S1 k_NL_PE = 100/h (free / not tabulated in source)
    k_ns_pe        <- fixed(280)         ; label("PLD forward rate k_NS_PE for NSPE (1/h)")                  # Table S1 k_NS_PE = 280/h; JBC 2004 279:5298-5305
    km_na_pe       <- fixed(2800)        ; label("PLD Km for NAPE Km_NA_PE (nM)")                            # Table S1 Km_NA_PE = 2800 nM; JBC 2004 279:5298-5305
    km_no_pe       <- fixed(2900)        ; label("PLD Km for NOPE Km_NO_PE (nM)")                            # Table S1 Km_NO_PE = 2900 nM
    km_np_pe       <- fixed(3300)        ; label("PLD Km for NPPE Km_NP_PE (nM)")                            # Table S1 Km_NP_PE = 3300 nM
    km_nl_pe       <- fixed(1000)        ; label("PLD Km for NLPE Km_NL_PE (nM; free/assumed)")              # Table S1 Km_NL_PE = 1000 nM (free)
    km_ns_pe       <- fixed(3400)        ; label("PLD Km for NSPE Km_NS_PE (nM)")                            # Table S1 Km_NS_PE = 3400 nM

    # ---------- PLD product-inhibition Ki (from ethanolamide products) --
    ki_a           <- fixed(230)         ; label("PLD Ki by AEA Ki_A (nM)")                                  # Table S1 Ki_A = 230 nM; Anal Biochem 2005 339:113-120 (mouse brain NAPE-PLD)
    ki_o           <- fixed(240)         ; label("PLD Ki by OEA Ki_O (nM)")                                  # Table S1 Ki_O = 240 nM
    ki_p           <- fixed(6700)        ; label("PLD Ki by PEA Ki_P (nM)")                                  # Table S1 Ki_P = 6700 nM
    ki_l           <- fixed(1000)        ; label("PLD Ki by LEA Ki_L (nM; free/assumed)")                    # Table S1 Ki_L = 1000 nM (free)
    ki_s           <- fixed(840)         ; label("PLD Ki by SEA Ki_S (nM)")                                  # Table S1 Ki_S = 840 nM

    # ---------- PLD concentrations in brain and ROB (SBML values) -------
    pld_brain      <- fixed(1e7)         ; label("PLD concentration in brain PLD_b (nM; SBML)")              # SBML Model S1 PLD_b = 1e7 nM (Errata: Table S1 lists 1e6 nM -- SBML value used because it produces the fitted CT profiles)
    pld_rob        <- fixed(1e7)         ; label("PLD concentration in ROB PLD_r (nM; SBML)")                # SBML Model S1 PLD_r = 1e7 nM (Errata: Table S1 lists 1e6 nM)

    # ---------- FAAH catalytic constant + relative substrate specificity
    kcat_faah      <- fixed(18000)       ; label("FAAH catalytic constant kcat_FAAH (1/h)")                  # Table S1 kcat_FAAH = 18000/h; Biochemistry 38:9804-9812 (rattus norvegicus brain)
    a_faah_a       <- fixed(1)           ; label("FAAH relative activity toward AEA a_FAAH_A (unitless)")    # Table S1 a_FAAH_A = 1 (reference)
    a_faah_o       <- fixed(5.7)         ; label("FAAH relative activity toward OEA a_FAAH_O (unitless)")    # Table S1 a_FAAH_O = 5.7; fitted Pfizer CT I
    a_faah_p       <- fixed(37.8)        ; label("FAAH relative activity toward PEA a_FAAH_P (unitless)")    # Table S1 a_FAAH_P = 37.8; fitted Pfizer CT I
    a_faah_l       <- fixed(1.15)        ; label("FAAH relative activity toward LEA a_FAAH_L (unitless)")    # Table S1 a_FAAH_L = 1.15; fitted Pfizer CT I
    a_faah_s       <- fixed(1)           ; label("FAAH relative activity toward SEA a_FAAH_S (unitless; free)") # Table S1 a_FAAH_S = 1 (free / assumed reference)

    # ---------- FAAH total pool + tissue distribution factors ----------
    faah_t         <- fixed(78)          ; label("FAAH total steady-state concentration FAAH_t (nM)")        # Table S1 FAAH_t = 78 nM; fitted Pfizer CT I
    b_faah_liver     <- fixed(1)         ; label("FAAH tissue factor liver b_FAAH_Liver (unitless)")         # Table S1 b_FAAH_Liver = 1 (reference)
    b_faah_brain     <- fixed(0.197)     ; label("FAAH tissue factor brain b_FAAH_Brain (unitless)")         # Table S1 b_FAAH_Brain = 0.197; calculated per BBA 1347:212-218
    b_faah_gut       <- fixed(0.034)     ; label("FAAH tissue factor gut b_FAAH_Gut (unitless)")             # Table S1 b_FAAH_Gut = 0.034
    b_faah_spleen    <- fixed(0.03)      ; label("FAAH tissue factor spleen b_FAAH_Spleen (unitless)")       # Table S1 b_FAAH_Spleen = 0.03
    b_faah_kidney    <- fixed(0.069)     ; label("FAAH tissue factor kidney b_FAAH_Kidney (unitless)")       # Table S1 b_FAAH_Kidney = 0.069
    b_faah_lungs     <- fixed(0.032)     ; label("FAAH tissue factor lungs b_FAAH_Lungs (unitless)")         # Table S1 b_FAAH_Lungs = 0.032
    b_faah_testis    <- fixed(0.126)     ; label("FAAH tissue factor testis b_FAAH_Testis (unitless)")       # Table S1 b_FAAH_Testis = 0.126
    b_faah_mec       <- fixed(0.137)     ; label("FAAH tissue factor MEC b_FAAH_MEC (unitless)")             # Table S1 b_FAAH_MEC = 0.137
    b_faah_leucocytes <- fixed(0)        ; label("FAAH tissue factor leucocytes b_FAAH_Leucocytes (unitless)") # SBML Model S1 b_FAAH_Leucocytes = 0

    # ---------- FAAH substrate Km values (rat brain in vitro) ----------
    km_faah_a      <- fixed(8200)        ; label("FAAH Km for AEA Km_FAAH_A (nM)")                           # Table S1 Km_FAAH_A = 8200 nM; JBC 270(11):6030-6035 (rat brain)
    km_faah_o      <- fixed(52200)       ; label("FAAH Km for OEA Km_FAAH_O (nM)")                           # Table S1 Km_FAAH_O = 52200 nM
    km_faah_p      <- fixed(543000)      ; label("FAAH Km for PEA Km_FAAH_P (nM)")                           # Table S1 Km_FAAH_P = 543000 nM
    km_faah_l      <- fixed(10800)       ; label("FAAH Km for LEA Km_FAAH_L (nM)")                           # Table S1 Km_FAAH_L = 10800 nM
    km_faah_s      <- fixed(10000)       ; label("FAAH Km for SEA Km_FAAH_S (nM; free/assumed)")             # Table S1 Km_FAAH_S = 10000 nM (free)

    # ---------- FAAH protein turnover + PF-04457845 irreversible inhib. -
    k_deg_faah     <- fixed(0.0051)      ; label("FAAH protein degradation rate k_deg_FAAH (1/h)")           # Table S1 k_deg_FAAH = 0.0051/h; fitted Pfizer CT I
    k_inh          <- fixed(1.1)         ; label("PF-04457845 irreversible-inhibition constant k_inh (1/(nM*h))") # Table S1 k_inh = 1.1 nM^-1/h; fitted Pfizer CT I

    # ---------- NAAA (FAAH-independent clearance) parameters ----------
    kcl_a          <- fixed(1.74)        ; label("NAAA clearance constant for AEA kcl_A (1/h)")              # Table S1 kcl_A = 1.74/h; fitted Pfizer CT I
    kcl_o          <- fixed(2.5)         ; label("NAAA clearance constant for OEA kcl_O (1/h)")              # Table S1 kcl_O = 2.5/h
    kcl_p          <- fixed(2.61)        ; label("NAAA clearance constant for PEA kcl_P (1/h)")              # Table S1 kcl_P = 2.61/h
    kcl_l          <- fixed(1.25)        ; label("NAAA clearance constant for LEA kcl_L (1/h)")              # Table S1 kcl_L = 1.25/h
    kcl_s          <- fixed(1.2)         ; label("NAAA clearance constant for SEA kcl_S (1/h; free)")        # Table S1 kcl_S = 1.2/h (free)

    # ---------- NAAA tissue distribution factors b_NAAA_<tissue> --------
    b_naaa_liver   <- fixed(1)           ; label("NAAA tissue factor liver (unitless)")                      # Table S1 b_NAAA_Liver = 1 (reference)
    b_naaa_gut     <- fixed(0.2)         ; label("NAAA tissue factor gut (unitless)")                        # Table S1 b_NAAA_Gut = 0.2; JBC 276(38):35552-35557
    b_naaa_spleen  <- fixed(8)           ; label("NAAA tissue factor spleen (unitless)")                     # Table S1 b_NAAA_Spleen = 8
    b_naaa_kidney  <- fixed(0.6)         ; label("NAAA tissue factor kidney (unitless)")                     # Table S1 b_NAAA_Kidney = 0.6
    b_naaa_heart   <- fixed(0.2)         ; label("NAAA tissue factor heart (unitless)")                      # Table S1 b_NAAA_Heart = 0.2
    b_naaa_lungs   <- fixed(14)          ; label("NAAA tissue factor lungs (unitless)")                      # Table S1 b_NAAA_Lungs = 14
    b_naaa_thymus  <- fixed(4)           ; label("NAAA tissue factor thymus (unitless)")                     # Table S1 b_NAAA_Thymus = 4
    b_naaa_testis  <- fixed(0.6)         ; label("NAAA tissue factor testis (unitless)")                     # Table S1 b_NAAA_Testis = 0.6
    # Table S1 b_NAAA_Brain = 0.6 is defined but never used in the SBML
    # rate laws: brain NAAA-mediated clearance in the SBML (vA_UE_b etc.)
    # uses b_FAAH_Brain as the tissue-activity scalar. Not encoded.

    # ---------- Ethanolamide transport parameters -----------------------
    # ROB <-> plasma: passive diffusion for OEA/PEA/LEA/SEA + AEA
    # transporter for AEA; MEC <-> plasma similarly; brain <-> MEC
    # analogous. Ktr_p_r_<X> = partition-coefficient-based reverse
    # transport constant calculated by Pfizer from logP; ktr_r_p is the
    # aggregate rate constant.
    ktr_r_p        <- fixed(100)         ; label("Ethanolamide transport rate ROB-plasma ktr_r_p (unitless factor)")  # SBML Model S1 ktr_r_p = 100 (Table S1 not explicitly listed)
    ktr_m_p_a      <- fixed(150)         ; label("AEA transporter rate MEC-plasma ktr_m_p_A (nM/h)")                  # Table S1 ktr_m_p_A = 150 nM/h; fitted per BBB model (Thromb Haemost 2006;95:117-127)
    ktr_m_p_o      <- fixed(10)          ; label("OEA transport rate MEC-plasma ktr_m_p_O (1/h; rapid-equilibrium)")  # Table S1 ktr_m_p_O = 10 (rapid-equilibrium fit)
    ktr_m_p_p      <- fixed(10)          ; label("PEA transport rate MEC-plasma ktr_m_p_P (1/h; rapid-equilibrium)")  # Table S1 ktr_m_p_P = 10
    ktr_m_p_l      <- fixed(0)           ; label("LEA transport rate MEC-plasma ktr_m_p_L (1/h; SBML)")               # SBML Model S1 ktr_m_p_L = 0 (Table S1 shows 10 fitted; SBML value used for consistency with published simulation)
    ktr_m_p_s      <- fixed(10)          ; label("SEA transport rate MEC-plasma ktr_m_p_S (1/h; rapid-equilibrium)")  # Table S1 ktr_m_p_S = 10
    km_p_m_a       <- fixed(1)           ; label("AEA transporter Km_p_m_A (nM)")                                     # Table S1 Km_p_m_A = 1 nM; fitted per BBB model
    ktr_p_m_a      <- fixed(1.89)        ; label("AEA membrane partition coefficient Ktr_p_m_A (unitless)")           # Table S1 Ktr_p_m_A = 1.89; calculated from logP by Pfizer
    ktr_p_m_o      <- fixed(9.07)        ; label("OEA membrane partition coefficient Ktr_p_m_O (unitless)")           # Table S1 Ktr_p_m_O = 9.07
    ktr_p_m_p      <- fixed(2.65)        ; label("PEA membrane partition coefficient Ktr_p_m_P (unitless)")           # Table S1 Ktr_p_m_P = 2.65
    ktr_p_m_l      <- fixed(2.77)        ; label("LEA membrane partition coefficient Ktr_p_m_L (unitless)")           # Table S1 Ktr_p_m_L = 2.77
    ktr_p_m_s      <- fixed(30.01)       ; label("SEA membrane partition coefficient Ktr_p_m_S (unitless)")           # Table S1 Ktr_p_m_S = 30.01
    ktr_p_r_a      <- fixed(0.62)        ; label("AEA plasma-ROB partition coefficient Ktr_p_r_A (unitless)")         # Table S1 Ktr_p_r_A = 0.62; calculated from logP
    ktr_p_r_o      <- fixed(2.8)         ; label("OEA plasma-ROB partition coefficient Ktr_p_r_O (unitless)")         # Table S1 Ktr_p_r_O = 2.8
    ktr_p_r_p      <- fixed(0.85)        ; label("PEA plasma-ROB partition coefficient Ktr_p_r_P (unitless)")         # Table S1 Ktr_p_r_P = 0.85
    ktr_p_r_l      <- fixed(0.89)        ; label("LEA plasma-ROB partition coefficient Ktr_p_r_L (unitless)")         # Table S1 Ktr_p_r_L = 0.89
    ktr_p_r_s      <- fixed(9.19)        ; label("SEA plasma-ROB partition coefficient Ktr_p_r_S (unitless)")         # Table S1 Ktr_p_r_S = 9.19

    # ---------- PF-04457845 PK parameters ------------------------------
    emax_pfm       <- fixed(0.773)       ; label("PF-04457845 bioavailability Emax Emax_PFM (fraction; SBML)")   # SBML Model S1 Emax_PFM = 0.773 (Errata: Table S1 lists 0.0773 -- SBML value used)
    ed50           <- fixed(0.53)        ; label("PF-04457845 bioavailability ED50 (mg)")                        # Table S1 ED50 = 0.53 mg; fitted Pfizer CT I
    kabs_pfm       <- fixed(2.2)         ; label("PF-04457845 absorption rate kabs_PFM (1/h)")                   # Table S1 kabs_PFM = 2.2/h; fitted Pfizer CT I
    kin_pfm        <- fixed(0.117)       ; label("PF-04457845 plasma-to-ROB rate kin_PFM (1/h)")                 # Table S1 kin_PFM = 0.117/h; fitted Pfizer CT I
    kout_pfm       <- fixed(0.18)        ; label("PF-04457845 ROB-to-plasma rate kout_PFM (1/h)")                # Table S1 kout_PFM = 0.18/h; fitted Pfizer CT I
    klinear_pfm    <- fixed(0.0803)      ; label("PF-04457845 linear elimination rate klinear_PFM (1/h)")        # Table S1 k_linear_PFM = 0.0803/h; fitted Pfizer CT I
    vm_pfm         <- fixed(1511)        ; label("PF-04457845 saturable elimination Vmax_PFM (ng/h)")            # Table S1 Vm_PFM = 1511 ng/h; fitted Pfizer CT I
    km_pfm         <- fixed(26.1)        ; label("PF-04457845 saturable elimination Km_PFM (ng/L)")              # Table S1 Km_PFM = 26.1 ng/L; fitted Pfizer CT I
    vss_pf         <- fixed(58.328)      ; label("PF-04457845 plasma distribution volume Vss_PFM (L; Errata)")   # Table S1 Vss_PFM = 58.328 (Errata: Table S1 label 'mL' is a supplement typo -- the model math and observed PK profile only balance with L)
    kp_b_pf        <- fixed(1.3)         ; label("PF-04457845 brain partition coefficient Kp_b_PF (unitless)")   # Table S1 Kp_b_PF = 1.3 (free)
    kp_r_pf        <- fixed(1.5)         ; label("PF-04457845 ROB partition coefficient Kp_r_PF (unitless; SBML)") # SBML Model S1 Kp_r_PF = 1.5 (Errata: Table S1 has a duplicate 'Kp_m_PF' row with the value 1.5 that is actually Kp_r_PF)
    kp_m_pf        <- fixed(1.3)         ; label("PF-04457845 MEC partition coefficient Kp_m_PF (unitless)")     # Table S1 Kp_m_PF = 1.3 (free)
    m_pf           <- fixed(455.4)       ; label("PF-04457845 molecular mass M_PF (g/mol)")                      # Table S1 M_PF = 455.4 g/mol; Pfizer data

    # ---------- Ethanolamide molecular masses (g/mol) -------------------
    m_a            <- fixed(347.5)       ; label("AEA molecular mass M_A (g/mol)")                               # Table S1 M_A = 347.5 g/mol
    m_o            <- fixed(325.5)       ; label("OEA molecular mass M_O (g/mol)")                               # Table S1 M_O = 325.5 g/mol
    m_p            <- fixed(299.5)       ; label("PEA molecular mass M_P (g/mol)")                               # Table S1 M_P = 299.5 g/mol
    m_l            <- fixed(323.5)       ; label("LEA molecular mass M_L (g/mol)")                               # Table S1 M_L = 323.5 g/mol
    m_s            <- fixed(321.5)       ; label("SEA molecular mass M_S (g/mol)")                               # Table S1 M_S = 321.5 g/mol

    # ---------- CB1 receptor binding parameters -------------------------
    kd_cb1_a       <- fixed(239.2)       ; label("CB1 in vitro Kd for AEA Kd_CB1_A (nM)")                        # Table S1 Kd_CB1_A = 239.2 nM; Br J Pharmacol 2007 152:583-593 (meta-analysis)

    # ---------- Residual error (FIXED at 0 -- paper reports no RUV) -----
    # Paper is a deterministic mechanistic model with no residual error
    # component; per operator policy for unreported RUV, encoded as
    # fixed(0). All ~ prop() / add() forms below use these zero SDs.
    propSd            <- fixed(0)        ; label("Proportional residual SD on PF-04457845 plasma Cc (not fitted)")
    propSd_Cc_AEA     <- fixed(0)        ; label("Proportional residual SD on AEA plasma Cc_AEA (not fitted)")
    propSd_Cc_OEA     <- fixed(0)        ; label("Proportional residual SD on OEA plasma Cc_OEA (not fitted)")
    propSd_Cc_PEA     <- fixed(0)        ; label("Proportional residual SD on PEA plasma Cc_PEA (not fitted)")
    propSd_Cc_LEA     <- fixed(0)        ; label("Proportional residual SD on LEA plasma Cc_LEA (not fitted)")
    propSd_Cc_SEA     <- fixed(0)        ; label("Proportional residual SD on SEA plasma Cc_SEA (not fitted)")
    addSd_FAAHact     <- fixed(0)        ; label("Additive residual SD on FAAHact (fraction; not fitted)")
    addSd_CB1occ      <- fixed(0)        ; label("Additive residual SD on CB1occ (fraction; not fitted)")
  })

  model({
    # =====================================================================
    # 1. PF-04457845 saturable oral bioavailability (Benson 2014
    #    supplement Table S1). The user provides dose_mg as a per-dose
    #    covariate and enters amt = dose_mg * 1e6 (ng) into pf_gut;
    #    F_PFM is applied via f(pf_gut) below.
    # =====================================================================
    f_pfm <- emax_pfm * DOSE / (ed50 + DOSE)

    # =====================================================================
    # 2. PF-04457845 concentrations by compartment (Benson 2014 supplement
    #    Table S1 assignment rules PF_p / PF_b / PF_r / PF_m). pf_p_amt is
    #    the drug amount in the plasma sub-pool of the drug distribution
    #    system (in ng); dividing by (M_PF * vss_pf) gives nmol/L (nM).
    #    The tissue concentrations follow the free-drug partitioning
    #    convention Kp_<tissue>_PF * plasma_free.
    # =====================================================================
    pf_p_conc  <- pf_p_amt / (m_pf * vss_pf)              # nM (nmol/L) in plasma sub-pool
    pf_b_conc  <- pf_p_conc * kp_b_pf                     # nM in brain
    pf_r_conc  <- pf_p_conc * kp_r_pf                     # nM in ROB
    pf_m_conc  <- pf_p_conc * kp_m_pf                     # nM in MEC

    # =====================================================================
    # 3. FAAH competitive-inhibition denominators (Benson 2014 supplement
    #    Table S1 FAAH_D_<compartment> assignment rules). Substrate
    #    competition among AEA/OEA/PEA/LEA/SEA at the FAAH active site.
    # =====================================================================
    faah_d_brain <- 1 + aea_brain/km_faah_a + oea_brain/km_faah_o + pea_brain/km_faah_p + lea_brain/km_faah_l + sea_brain/km_faah_s
    faah_d_rob   <- 1 + aea_rob/km_faah_a   + oea_rob/km_faah_o   + pea_rob/km_faah_p   + lea_rob/km_faah_l   + sea_rob/km_faah_s
    faah_d_mec   <- 1 + aea_mec/km_faah_a   + oea_mec/km_faah_o   + pea_mec/km_faah_p   + lea_mec/km_faah_l   + sea_mec/km_faah_s

    # =====================================================================
    # 4. NAPE-PLD competitive-inhibition denominators. NAPE and other
    #    N-acyl-PE precursors compete for the same PLD active site
    #    (slag1_<c>); ethanolamide products competitively inhibit the
    #    reverse (slag2_<c>). den_<c> = 1 + slag1_<c> + slag2_<c>.
    # =====================================================================
    slag1_brain <- nape_brain/km_na_pe + nope_brain/km_no_pe + nppe_brain/km_np_pe + nlpe_brain/km_nl_pe + nspe_brain/km_ns_pe
    slag2_brain <- aea_brain/ki_a + oea_brain/ki_o + pea_brain/ki_p + lea_brain/ki_l + sea_brain/ki_s
    den_brain   <- 1 + slag1_brain + slag2_brain

    slag1_rob   <- nape_rob/km_na_pe + nope_rob/km_no_pe + nppe_rob/km_np_pe + nlpe_rob/km_nl_pe + nspe_rob/km_ns_pe
    slag2_rob   <- aea_rob/ki_a + oea_rob/ki_o + pea_rob/ki_p + lea_rob/ki_l + sea_rob/ki_s
    den_rob     <- 1 + slag1_rob + slag2_rob

    # =====================================================================
    # 5. Composite tissue-weighted NAT / FAAH / NAAA activities acting on
    #    the ROB pool (Benson 2014 supplement Table S1 c_NAT_ROB,
    #    c_FAAH_ROB, c_NAAA_ROB assignment rules; sum over tissues of
    #    tissue_volume * relative_tissue_factor).
    # =====================================================================
    c_nat_rob  <- pancreas_vol*b_nat_pancreas + kidney_vol*b_nat_kidney + heart_vol*b_nat_heart + lungs_vol*b_nat_lungs + muscles_vol*b_nat_muscles + testis_vol*b_nat_testis + leucocytes_vol*b_nat_leucocytes
    c_faah_rob <- liver_vol*b_faah_liver + gut_vol*b_faah_gut + spleen_vol*b_faah_spleen + kidney_vol*b_faah_kidney + lungs_vol*b_faah_lungs + testis_vol*b_faah_testis + leucocytes_vol*b_faah_leucocytes
    c_naaa_rob <- liver_vol*b_naaa_liver + gut_vol*b_naaa_gut + spleen_vol*b_naaa_spleen + kidney_vol*b_naaa_kidney + heart_vol*b_naaa_heart + lungs_vol*b_naaa_lungs + thymus_vol*b_naaa_thymus + testis_vol*b_naaa_testis

    # =====================================================================
    # 6. Reaction rate expressions (Benson 2014 supplement Table S1
    #    columns "Rate law"). Each expression is in amount/time units
    #    (nmol/h for ethanolamides / NAPEs / FAAH; ng/h for the drug).
    #    ODE assembly below divides each amount-rate by the receiving
    #    compartment volume to update the concentration state.
    # =====================================================================

    # ---- BRAIN: FAAH-mediated ethanolamide degradation
    v_a_degr_brain <- brain_vol * faah_brain * kcat_faah * a_faah_a * aea_brain / (km_faah_a * faah_d_brain)
    v_o_degr_brain <- brain_vol * faah_brain * kcat_faah * a_faah_o * oea_brain / (km_faah_o * faah_d_brain)
    v_p_degr_brain <- brain_vol * faah_brain * kcat_faah * a_faah_p * pea_brain / (km_faah_p * faah_d_brain)
    v_l_degr_brain <- brain_vol * faah_brain * kcat_faah * a_faah_l * lea_brain / (km_faah_l * faah_d_brain)
    v_s_degr_brain <- brain_vol * faah_brain * kcat_faah * a_faah_s * sea_brain / (km_faah_s * faah_d_brain)

    # ---- BRAIN: NAPE synthesis by N-acyltransferase
    v_nape_syn_brain <- brain_vol * vmax_nat * p_a * a_nat_a * b_nat_brain
    v_nope_syn_brain <- brain_vol * vmax_nat * p_o * a_nat_o * b_nat_brain
    v_nppe_syn_brain <- brain_vol * vmax_nat * p_p * a_nat_p * b_nat_brain
    v_nlpe_syn_brain <- brain_vol * vmax_nat * p_l * a_nat_l * b_nat_brain
    v_nspe_syn_brain <- brain_vol * vmax_nat * p_s * a_nat_s * b_nat_brain

    # ---- BRAIN: PLD hydrolysis of NAPEs into ethanolamides
    v_a_syn_brain <- brain_vol * pld_brain * k_na_pe * nape_brain / km_na_pe / den_brain
    v_o_syn_brain <- brain_vol * pld_brain * k_no_pe * nope_brain / km_no_pe / den_brain
    v_p_syn_brain <- brain_vol * pld_brain * k_np_pe * nppe_brain / km_np_pe / den_brain
    v_l_syn_brain <- brain_vol * pld_brain * k_nl_pe * nlpe_brain / km_nl_pe / den_brain
    v_s_syn_brain <- brain_vol * pld_brain * k_ns_pe * nspe_brain / km_ns_pe / den_brain

    # ---- BRAIN: FAAH protein turnover and PF-04457845 irreversible inhibition
    v_faah_syn_brain      <- brain_vol * faah_t * b_faah_brain * k_deg_faah
    v_faah_degr_brain     <- brain_vol * k_deg_faah * faah_brain
    v_faah_inh_brain      <- brain_vol * k_inh * faah_brain * pf_b_conc
    v_faah_inh_degr_brain <- brain_vol * k_deg_faah * faahinh_brain

    # ---- BRAIN: NAAA-mediated (unknown-enzyme) ethanolamide clearance
    v_a_ue_brain <- brain_vol * b_faah_brain * kcl_a * aea_brain
    v_o_ue_brain <- brain_vol * b_faah_brain * kcl_o * oea_brain
    v_p_ue_brain <- brain_vol * b_faah_brain * kcl_p * pea_brain
    v_l_ue_brain <- brain_vol * b_faah_brain * kcl_l * lea_brain
    v_s_ue_brain <- brain_vol * b_faah_brain * kcl_s * sea_brain

    # ---- ROB: FAAH-mediated degradation
    v_a_degr_rob <- rob_vol * faah_rob * kcat_faah * a_faah_a * aea_rob / (km_faah_a * faah_d_rob)
    v_o_degr_rob <- rob_vol * faah_rob * kcat_faah * a_faah_o * oea_rob / (km_faah_o * faah_d_rob)
    v_p_degr_rob <- rob_vol * faah_rob * kcat_faah * a_faah_p * pea_rob / (km_faah_p * faah_d_rob)
    v_l_degr_rob <- rob_vol * faah_rob * kcat_faah * a_faah_l * lea_rob / (km_faah_l * faah_d_rob)
    v_s_degr_rob <- rob_vol * faah_rob * kcat_faah * a_faah_s * sea_rob / (km_faah_s * faah_d_rob)

    # ---- ROB: NAPE synthesis (uses composite c_nat_rob rather than a single tissue factor)
    v_nape_syn_rob <- vmax_nat * p_a * a_nat_a * c_nat_rob
    v_nope_syn_rob <- vmax_nat * p_o * a_nat_o * c_nat_rob
    v_nppe_syn_rob <- vmax_nat * p_p * a_nat_p * c_nat_rob
    v_nlpe_syn_rob <- vmax_nat * p_l * a_nat_l * c_nat_rob
    v_nspe_syn_rob <- vmax_nat * p_s * a_nat_s * c_nat_rob

    # ---- ROB: PLD hydrolysis
    v_a_syn_rob <- rob_vol * pld_rob * k_na_pe * nape_rob / km_na_pe / den_rob
    v_o_syn_rob <- rob_vol * pld_rob * k_no_pe * nope_rob / km_no_pe / den_rob
    v_p_syn_rob <- rob_vol * pld_rob * k_np_pe * nppe_rob / km_np_pe / den_rob
    v_l_syn_rob <- rob_vol * pld_rob * k_nl_pe * nlpe_rob / km_nl_pe / den_rob
    v_s_syn_rob <- rob_vol * pld_rob * k_ns_pe * nspe_rob / km_ns_pe / den_rob

    # ---- ROB: FAAH turnover and inhibition
    v_faah_syn_rob      <- faah_t * c_faah_rob * k_deg_faah
    v_faah_degr_rob     <- rob_vol * k_deg_faah * faah_rob
    v_faah_inh_rob      <- rob_vol * k_inh * faah_rob * pf_r_conc
    v_faah_inh_degr_rob <- rob_vol * k_deg_faah * faahinh_rob

    # ---- ROB: NAAA clearance
    v_a_ue_rob <- c_naaa_rob * kcl_a * aea_rob
    v_o_ue_rob <- c_naaa_rob * kcl_o * oea_rob
    v_p_ue_rob <- c_naaa_rob * kcl_p * pea_rob
    v_l_ue_rob <- c_naaa_rob * kcl_l * lea_rob
    v_s_ue_rob <- c_naaa_rob * kcl_s * sea_rob

    # ---- MEC: FAAH-mediated degradation
    v_a_degr_mec <- mec_vol * faah_mec * kcat_faah * a_faah_a * aea_mec / (km_faah_a * faah_d_mec)
    v_o_degr_mec <- mec_vol * faah_mec * kcat_faah * a_faah_o * oea_mec / (km_faah_o * faah_d_mec)
    v_p_degr_mec <- mec_vol * faah_mec * kcat_faah * a_faah_p * pea_mec / (km_faah_p * faah_d_mec)
    v_l_degr_mec <- mec_vol * faah_mec * kcat_faah * a_faah_l * lea_mec / (km_faah_l * faah_d_mec)
    v_s_degr_mec <- mec_vol * faah_mec * kcat_faah * a_faah_s * sea_mec / (km_faah_s * faah_d_mec)

    # ---- MEC: FAAH turnover and inhibition (no NAPE synthesis; no NAAA clearance in MEC per SBML)
    v_faah_syn_mec      <- mec_vol * faah_t * b_faah_mec * k_deg_faah
    v_faah_degr_mec     <- mec_vol * k_deg_faah * faah_mec
    v_faah_inh_mec      <- mec_vol * k_inh * faah_mec * pf_m_conc
    v_faah_inh_degr_mec <- mec_vol * k_deg_faah * faahinh_mec

    # ---- Ethanolamide transport (ROB -> plasma, PLASMA multiplier)
    v_a_r_p <- plasma_vol * ktr_r_p * (aea_rob - aea_plasma * ktr_p_r_a) / (aea_rob + aea_plasma + km_p_m_a)
    v_o_r_p <- plasma_vol * ktr_r_p * (oea_rob - oea_plasma * ktr_p_r_o)
    v_p_r_p <- plasma_vol * ktr_r_p * (pea_rob - pea_plasma * ktr_p_r_p)
    v_l_r_p <- plasma_vol * ktr_r_p * (lea_rob - lea_plasma * ktr_p_r_l)
    v_s_r_p <- plasma_vol * ktr_r_p * (sea_rob - sea_plasma * ktr_p_r_s)

    # ---- Ethanolamide transport (MEC -> plasma, MEC multiplier)
    v_a_m_p <- mec_vol * ktr_m_p_a * (aea_mec - aea_plasma * ktr_p_m_a) / (aea_mec + aea_plasma + km_p_m_a)
    v_o_m_p <- mec_vol * ktr_m_p_o * (oea_mec - oea_plasma * ktr_p_m_o)
    v_p_m_p <- mec_vol * ktr_m_p_p * (pea_mec - pea_plasma * ktr_p_m_p)
    v_l_m_p <- mec_vol * ktr_m_p_l * (lea_mec - lea_plasma * ktr_p_m_l)
    v_s_m_p <- mec_vol * ktr_m_p_s * (sea_mec - sea_plasma * ktr_p_m_s)

    # ---- Ethanolamide transport (BRAIN -> MEC, MEC multiplier)
    v_a_b_m <- mec_vol * ktr_m_p_a * (aea_brain - aea_mec) / (aea_mec + aea_brain + km_p_m_a)
    v_o_b_m <- mec_vol * ktr_m_p_o * (oea_brain - oea_mec)
    v_p_b_m <- mec_vol * ktr_m_p_p * (pea_brain - pea_mec)
    v_l_b_m <- mec_vol * ktr_m_p_l * (lea_brain - lea_mec)
    v_s_b_m <- mec_vol * ktr_m_p_s * (sea_brain - sea_mec)

    # ---- PF-04457845 disposition (in ng amount units; Default compartment size 1 L)
    v_absorp <- kabs_pfm * pf_gut
    v_dist   <- kout_pfm * pf_p_amt - kin_pfm * pf_r_amt
    v_elim   <- klinear_pfm * pf_p_amt + vm_pfm * pf_p_amt / (km_pfm + pf_p_amt / vss_pf) / vss_pf

    # =====================================================================
    # 7. ODE system. States are CONCENTRATIONS (nM) for ethanolamides /
    #    NAPEs / FAAH / FAAHinh; concentration change = (amount-rate) /
    #    compartment_volume. PF-04457845 states pf_gut / pf_p_amt /
    #    pf_r_amt are AMOUNTS in ng (Default compartment size = 1 L).
    # =====================================================================

    # -- BRAIN (12) ---------------------------------------------------
    d/dt(aea_brain)  <- (v_a_syn_brain - v_a_degr_brain - v_a_ue_brain - v_a_b_m) / brain_vol
    d/dt(oea_brain)  <- (v_o_syn_brain - v_o_degr_brain - v_o_ue_brain - v_o_b_m) / brain_vol
    d/dt(pea_brain)  <- (v_p_syn_brain - v_p_degr_brain - v_p_ue_brain - v_p_b_m) / brain_vol
    d/dt(lea_brain)  <- (v_l_syn_brain - v_l_degr_brain - v_l_ue_brain - v_l_b_m) / brain_vol
    d/dt(sea_brain)  <- (v_s_syn_brain - v_s_degr_brain - v_s_ue_brain - v_s_b_m) / brain_vol
    d/dt(nape_brain) <- (v_nape_syn_brain - v_a_syn_brain) / brain_vol
    d/dt(nope_brain) <- (v_nope_syn_brain - v_o_syn_brain) / brain_vol
    d/dt(nppe_brain) <- (v_nppe_syn_brain - v_p_syn_brain) / brain_vol
    d/dt(nlpe_brain) <- (v_nlpe_syn_brain - v_l_syn_brain) / brain_vol
    d/dt(nspe_brain) <- (v_nspe_syn_brain - v_s_syn_brain) / brain_vol
    d/dt(faah_brain) <- (v_faah_syn_brain - v_faah_degr_brain - v_faah_inh_brain) / brain_vol
    d/dt(faahinh_brain) <- (v_faah_inh_brain - v_faah_inh_degr_brain) / brain_vol

    # -- ROB (12) -----------------------------------------------------
    d/dt(aea_rob)  <- (v_a_syn_rob - v_a_degr_rob - v_a_ue_rob - v_a_r_p) / rob_vol
    d/dt(oea_rob)  <- (v_o_syn_rob - v_o_degr_rob - v_o_ue_rob - v_o_r_p) / rob_vol
    d/dt(pea_rob)  <- (v_p_syn_rob - v_p_degr_rob - v_p_ue_rob - v_p_r_p) / rob_vol
    d/dt(lea_rob)  <- (v_l_syn_rob - v_l_degr_rob - v_l_ue_rob - v_l_r_p) / rob_vol
    d/dt(sea_rob)  <- (v_s_syn_rob - v_s_degr_rob - v_s_ue_rob - v_s_r_p) / rob_vol
    d/dt(nape_rob) <- (v_nape_syn_rob - v_a_syn_rob) / rob_vol
    d/dt(nope_rob) <- (v_nope_syn_rob - v_o_syn_rob) / rob_vol
    d/dt(nppe_rob) <- (v_nppe_syn_rob - v_p_syn_rob) / rob_vol
    d/dt(nlpe_rob) <- (v_nlpe_syn_rob - v_l_syn_rob) / rob_vol
    d/dt(nspe_rob) <- (v_nspe_syn_rob - v_s_syn_rob) / rob_vol
    d/dt(faah_rob) <- (v_faah_syn_rob - v_faah_degr_rob - v_faah_inh_rob) / rob_vol
    d/dt(faahinh_rob) <- (v_faah_inh_rob - v_faah_inh_degr_rob) / rob_vol

    # -- MEC (7) ------------------------------------------------------
    d/dt(aea_mec) <- (v_a_b_m - v_a_degr_mec - v_a_m_p) / mec_vol
    d/dt(oea_mec) <- (v_o_b_m - v_o_degr_mec - v_o_m_p) / mec_vol
    d/dt(pea_mec) <- (v_p_b_m - v_p_degr_mec - v_p_m_p) / mec_vol
    d/dt(lea_mec) <- (v_l_b_m - v_l_degr_mec - v_l_m_p) / mec_vol
    d/dt(sea_mec) <- (v_s_b_m - v_s_degr_mec - v_s_m_p) / mec_vol
    d/dt(faah_mec) <- (v_faah_syn_mec - v_faah_degr_mec - v_faah_inh_mec) / mec_vol
    d/dt(faahinh_mec) <- (v_faah_inh_mec - v_faah_inh_degr_mec) / mec_vol

    # -- PLASMA (5) --------------------------------------------------
    d/dt(aea_plasma) <- (v_a_m_p + v_a_r_p) / plasma_vol
    d/dt(oea_plasma) <- (v_o_m_p + v_o_r_p) / plasma_vol
    d/dt(pea_plasma) <- (v_p_m_p + v_p_r_p) / plasma_vol
    d/dt(lea_plasma) <- (v_l_m_p + v_l_r_p) / plasma_vol
    d/dt(sea_plasma) <- (v_s_m_p + v_s_r_p) / plasma_vol

    # -- PF-04457845 drug states (3; amounts in ng) ------------------
    d/dt(pf_gut)   <- -v_absorp
    d/dt(pf_p_amt) <-  v_absorp - v_dist - v_elim
    d/dt(pf_r_amt) <-  v_dist

    # =====================================================================
    # 8. Initial conditions (Benson 2014 supplement Model S1 SBML export;
    #    reproduces the steady-state biomarker levels at t=0). All
    #    ethanolamide, NAPE, FAAH and FAAHinh initial concentrations are
    #    fixed at the SBML-reported steady-state values; drug states
    #    start at 0.
    # =====================================================================
    # -- BRAIN
    aea_brain(0)   <- 0.7493309
    oea_brain(0)   <- 20.77858
    pea_brain(0)   <- 6.541209
    lea_brain(0)   <- 2.319571
    sea_brain(0)   <- 3.427807
    nape_brain(0)  <- 3.879041e-5
    nope_brain(0)  <- 8.814287e-4
    nppe_brain(0)  <- 1.732296e-4
    nlpe_brain(0)  <- 7.550331e-5
    nspe_brain(0)  <- 1.272629e-4
    faah_brain(0)  <- 15.366
    faahinh_brain(0) <- 0
    # -- ROB
    aea_rob(0)     <- 0.5419204
    oea_rob(0)     <- 14.23822
    pea_rob(0)     <- 4.121915
    lea_rob(0)     <- 1.705466
    sea_rob(0)     <- 2.515968
    nape_rob(0)    <- 4.241633e-6
    nope_rob(0)    <- 9.638198e-5
    nppe_rob(0)    <- 1.894222e-5
    nlpe_rob(0)    <- 8.256095e-6
    nspe_rob(0)    <- 1.391587e-5
    faah_rob(0)    <- 2.165868
    faahinh_rob(0) <- 0
    # -- MEC
    aea_mec(0)     <- 0.97761
    oea_mec(0)     <- 16.3219
    pea_mec(0)     <- 5.809415
    lea_mec(0)     <- 0
    sea_mec(0)     <- 2.968774
    faah_mec(0)    <- 10.686
    faahinh_mec(0) <- 0
    # -- PLASMA
    aea_plasma(0)  <- 0.8740574
    oea_plasma(0)  <- 5.085073
    pea_plasma(0)  <- 4.849307
    lea_plasma(0)  <- 1.916254
    sea_plasma(0)  <- 0.273772
    # -- Drug states start at zero
    pf_gut(0)      <- 0
    pf_p_amt(0)    <- 0
    pf_r_amt(0)    <- 0

    # =====================================================================
    # 9. Observations (matching paper Figures 3-5).
    #    * Cc      = PF-04457845 plasma concentration in ng/mL (SBML PFG_p)
    #    * Cc_AEA  / Cc_OEA / Cc_PEA / Cc_LEA / Cc_SEA = ethanolamide
    #      plasma concentrations in ng/mL (SBML AG_p / OG_p / PG_p / LG_p /
    #      SG_p; converted from plasma nM via molecular mass)
    #    * FAAHact = ROB-pool FAAH activity as fraction of baseline (0-1);
    #      captures the paper's plasma-leukocyte FAAH assay in Figure 4a,b
    #    * CB1occ  = brain CB1 receptor occupancy fraction (0-1); the
    #      paper's Figure 5. 2-AG contribution drops out because
    #      ag2_brain is fixed at 0 per paper scope.
    # =====================================================================
    Cc       <- 1e-3 * pf_p_amt / vss_pf                                       # ng/mL
    Cc_AEA   <- 1e-3 * aea_plasma * m_a                                        # ng/mL from nM * g/mol
    Cc_OEA   <- 1e-3 * oea_plasma * m_o
    Cc_PEA   <- 1e-3 * pea_plasma * m_p
    Cc_LEA   <- 1e-3 * lea_plasma * m_l
    Cc_SEA   <- 1e-3 * sea_plasma * m_s
    FAAHact  <- faah_rob / (faah_rob + faahinh_rob)                            # fraction
    CB1occ   <- (aea_brain / kd_cb1_a) / (1 + aea_brain / kd_cb1_a)            # fraction

    # =====================================================================
    # 10. Bioavailability of the oral PF-04457845 dose (saturable Emax
    #     from Benson 2014 supplement Table S1). The user is expected to
    #     supply amt = dose_mg * 1e6 (ng) when dosing pf_gut; F_PFM
    #     reduces the entering amount according to dose_mg.
    # =====================================================================
    f(pf_gut) <- f_pfm

    # =====================================================================
    # 11. Residual-error models. All SDs are fixed at 0 because the
    #     paper reports no residual variability (deterministic
    #     mechanistic model). Users who wish to add plausible RUV for
    #     Monte-Carlo simulation should override the propSd_* and
    #     addSd_* parameters through rxode2's model-piping helpers.
    # =====================================================================
    Cc       ~ prop(propSd)
    Cc_AEA   ~ prop(propSd_Cc_AEA)
    Cc_OEA   ~ prop(propSd_Cc_OEA)
    Cc_PEA   ~ prop(propSd_Cc_PEA)
    Cc_LEA   ~ prop(propSd_Cc_LEA)
    Cc_SEA   ~ prop(propSd_Cc_SEA)
    FAAHact  ~ add(addSd_FAAHact)
    CB1occ   ~ add(addSd_CB1occ)
  })
}
