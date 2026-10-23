Yang_2021_metformin_ddi_pbpk <- function() {
  description <- paste(
    "PBPK (whole-body, transporter-mediated drug-drug interaction;",
    "Phoenix WinNonlin 8.1). The Yang 2021 metformin model solved",
    "simultaneously with one co-administered perpetrator that competitively",
    "inhibits the metformin carriers: intestinal OCT3 (apical enterocyte),",
    "hepatic OCT1 (sinusoidal), renal OCT2 (tubule basolateral) and renal",
    "MATE1/2-K (tubule apical). Metformin absorption is nonlinear:",
    "Michaelis-Menten uptake by intestinal OCT3 and PMAT in duodenum,",
    "jejunum and ileum (regional expression and pH scaling) plus 30%",
    "passive, over a six-segment gut-lumen transit chain. The liver is",
    "split into hepatic blood and hepatocytes (OCT1 uptake, passive",
    "diffusion, metabolism); the kidney into renal blood and tubule cells",
    "(OCT2 uptake with back-flux, MATE efflux, glomerular filtration).",
    "The perpetrator slot is parameterised: cimetidine uses the same",
    "transporter-resolved liver/kidney structure (tr = 1), while the other",
    "five perpetrators (pyrimethamine, trimethoprim, ondansetron,",
    "rabeprazole, verapamil) use well-stirred liver and kidney (tr = 0).",
    "Dosing the perpetrator states with zero amount recovers the",
    "metformin-alone arm, so one model object gives both arms of each",
    "interaction. The model ships set to cimetidine; the drugLibrary",
    "metadata holds all six perpetrators' inputs and the validation",
    "vignette swaps them in to reproduce Table 2, Table 3, Table S4 and",
    "Figure 6. Deterministic typical-value model: the paper's VPC",
    "estimated inter-individual variances on Vmax,OCT3, CLint,up,OCT1,",
    "CLint,up,OCT2, CLint,eff,MATE and Kti but did not report them, so no",
    "IIV is implemented. Several values are as-run reproductions of the",
    "WinNonlin model that deviate from the printed equations/constants",
    "(the metformin M-M absorption constants, the back-solved hepatic",
    "passive-diffusion clearance, and the rabeprazole metabolic",
    "clearance); see the vignette Errata.",
    sep = " "
  )
  reference <- paste(
    "Yang Y, Zhang Z, Li P, Kong W, Liu X, Liu L. A Whole-Body",
    "Physiologically Based Pharmacokinetic Model Characterizing Interplay",
    "of OCTs and MATEs in Intestine, Liver and Kidney to Predict",
    "Drug-Drug Interactions of Metformin with Perpetrators. Pharmaceutics.",
    "2021;13(5):698. doi:10.3390/pharmaceutics13050698. Equations 1-17,",
    "Table 1 (drug parameters), Supplementary Table S1 (physiology),",
    "Table S2 (tissue:plasma ratios), Table S4 (metformin dose-ranging",
    "predictions) and Table S6 (DDI schedules).",
    sep = " "
  )
  vignette <- "Yang_2021_metformin_ddi_pbpk"
  units <- list(time = "min", dosing = "ug", concentration = "ug/mL")

  # Non-canonical states: the six-segment gut lumen and gut wall, the
  # liver blood/hepatocyte and kidney blood/tubule splits, and every state
  # of the perpetrator slot (bare organ + _perpetrator suffix, the
  # parameterised interacting-drug role from compartment-names.md).
  paper_specific_compartments <- c(
    "gw_stomach",
    "gw_duodenum",
    "gw_jejunum",
    "gw_ileum",
    "gw_cecum",
    "gw_colon",
    "liver_blood",
    "hepatocyte",
    "kidney_blood",
    "tubule",
    "stomach_perpetrator",
    "duodenum_perpetrator",
    "jejunum_perpetrator",
    "ileum_perpetrator",
    "cecum_perpetrator",
    "colon_perpetrator",
    "gw_stomach_perpetrator",
    "gw_duodenum_perpetrator",
    "gw_jejunum_perpetrator",
    "gw_ileum_perpetrator",
    "gw_cecum_perpetrator",
    "gw_colon_perpetrator",
    "spleen_perpetrator",
    "liver_blood_perpetrator",
    "hepatocyte_perpetrator",
    "liver_perpetrator",
    "kidney_blood_perpetrator",
    "tubule_perpetrator",
    "kidney_perpetrator",
    "lung_perpetrator",
    "arterial_perpetrator",
    "venous_perpetrator",
    "heart_perpetrator",
    "brain_perpetrator",
    "muscle_perpetrator",
    "adipose_perpetrator",
    "skin_perpetrator",
    "other_perpetrator"
  )

  compartmentData <- list(
    stomach = list(analyte = "metformin", units = "ug", specimen = "administration site", verified = TRUE),
    duodenum = list(analyte = "metformin", units = "ug", specimen = "administration site", verified = TRUE),
    jejunum = list(analyte = "metformin", units = "ug", specimen = "administration site", verified = TRUE),
    ileum = list(analyte = "metformin", units = "ug", specimen = "administration site", verified = TRUE),
    cecum = list(analyte = "metformin", units = "ug", specimen = "administration site", verified = TRUE),
    colon = list(analyte = "metformin", units = "ug", specimen = "administration site", verified = TRUE),
    a_feces = list(analyte = "metformin", units = "ug", specimen = "faeces", verified = TRUE),
    gw_stomach = list(analyte = "metformin", units = "ug", specimen = "tissue", verified = TRUE),
    gw_duodenum = list(analyte = "metformin", units = "ug", specimen = "tissue", verified = TRUE),
    gw_jejunum = list(analyte = "metformin", units = "ug", specimen = "tissue", verified = TRUE),
    gw_ileum = list(analyte = "metformin", units = "ug", specimen = "tissue", verified = TRUE),
    gw_cecum = list(analyte = "metformin", units = "ug", specimen = "tissue", verified = TRUE),
    gw_colon = list(analyte = "metformin", units = "ug", specimen = "tissue", verified = TRUE),
    spleen = list(analyte = "metformin", units = "ug", specimen = "tissue", verified = TRUE),
    liver_blood = list(analyte = "metformin", units = "ug", specimen = "whole blood", verified = TRUE),
    hepatocyte = list(analyte = "metformin", units = "ug", specimen = "tissue", verified = TRUE),
    kidney_blood = list(analyte = "metformin", units = "ug", specimen = "whole blood", verified = TRUE),
    tubule = list(analyte = "metformin", units = "ug", specimen = "tissue", verified = TRUE),
    a_urine = list(analyte = "metformin", units = "ug", specimen = "urine", verified = TRUE),
    lung = list(analyte = "metformin", units = "ug", specimen = "tissue", verified = TRUE),
    arterial = list(analyte = "metformin", units = "ug", specimen = "whole blood", verified = TRUE),
    venous = list(analyte = "metformin", units = "ug", specimen = "whole blood", verified = TRUE),
    heart = list(analyte = "metformin", units = "ug", specimen = "tissue", verified = TRUE),
    brain = list(analyte = "metformin", units = "ug", specimen = "tissue", verified = TRUE),
    muscle = list(analyte = "metformin", units = "ug", specimen = "tissue", verified = TRUE),
    adipose = list(analyte = "metformin", units = "ug", specimen = "tissue", verified = TRUE),
    skin = list(analyte = "metformin", units = "ug", specimen = "tissue", verified = TRUE),
    other = list(analyte = "metformin", units = "ug", specimen = "tissue", verified = TRUE),
    stomach_perpetrator = list(analyte = "cimetidine", units = "ug", specimen = "administration site", verified = TRUE),
    duodenum_perpetrator = list(
      analyte = "cimetidine",
      units = "ug",
      specimen = "administration site",
      verified = TRUE
    ),
    jejunum_perpetrator = list(analyte = "cimetidine", units = "ug", specimen = "administration site", verified = TRUE),
    ileum_perpetrator = list(analyte = "cimetidine", units = "ug", specimen = "administration site", verified = TRUE),
    cecum_perpetrator = list(analyte = "cimetidine", units = "ug", specimen = "administration site", verified = TRUE),
    colon_perpetrator = list(analyte = "cimetidine", units = "ug", specimen = "administration site", verified = TRUE),
    gw_stomach_perpetrator = list(analyte = "cimetidine", units = "ug", specimen = "tissue", verified = TRUE),
    gw_duodenum_perpetrator = list(analyte = "cimetidine", units = "ug", specimen = "tissue", verified = TRUE),
    gw_jejunum_perpetrator = list(analyte = "cimetidine", units = "ug", specimen = "tissue", verified = TRUE),
    gw_ileum_perpetrator = list(analyte = "cimetidine", units = "ug", specimen = "tissue", verified = TRUE),
    gw_cecum_perpetrator = list(analyte = "cimetidine", units = "ug", specimen = "tissue", verified = TRUE),
    gw_colon_perpetrator = list(analyte = "cimetidine", units = "ug", specimen = "tissue", verified = TRUE),
    spleen_perpetrator = list(analyte = "cimetidine", units = "ug", specimen = "tissue", verified = TRUE),
    liver_blood_perpetrator = list(analyte = "cimetidine", units = "ug", specimen = "whole blood", verified = TRUE),
    hepatocyte_perpetrator = list(analyte = "cimetidine", units = "ug", specimen = "tissue", verified = TRUE),
    liver_perpetrator = list(analyte = "cimetidine", units = "ug", specimen = "tissue", verified = TRUE),
    kidney_blood_perpetrator = list(analyte = "cimetidine", units = "ug", specimen = "whole blood", verified = TRUE),
    tubule_perpetrator = list(analyte = "cimetidine", units = "ug", specimen = "tissue", verified = TRUE),
    kidney_perpetrator = list(analyte = "cimetidine", units = "ug", specimen = "tissue", verified = TRUE),
    lung_perpetrator = list(analyte = "cimetidine", units = "ug", specimen = "tissue", verified = TRUE),
    arterial_perpetrator = list(analyte = "cimetidine", units = "ug", specimen = "whole blood", verified = TRUE),
    venous_perpetrator = list(analyte = "cimetidine", units = "ug", specimen = "whole blood", verified = TRUE),
    heart_perpetrator = list(analyte = "cimetidine", units = "ug", specimen = "tissue", verified = TRUE),
    brain_perpetrator = list(analyte = "cimetidine", units = "ug", specimen = "tissue", verified = TRUE),
    muscle_perpetrator = list(analyte = "cimetidine", units = "ug", specimen = "tissue", verified = TRUE),
    adipose_perpetrator = list(analyte = "cimetidine", units = "ug", specimen = "tissue", verified = TRUE),
    skin_perpetrator = list(analyte = "cimetidine", units = "ug", specimen = "tissue", verified = TRUE),
    other_perpetrator = list(analyte = "cimetidine", units = "ug", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species = "human",
    n_subjects = NA_integer_,
    n_studies = NA_integer_,
    age_range = "healthy adults (young, healthy volunteers in the cited DDI studies)",
    weight_range = "70 kg reference adult used for the physiology (Table S1)",
    disease_state = "healthy volunteers",
    dose_range = paste(
      "Oral metformin 0.1 mg microdose and 250-1000 mg; perpetrators oral",
      "cimetidine 300-800 mg, pyrimethamine 50-75 mg, trimethoprim",
      "200-1400 mg, ondansetron 8 mg, rabeprazole 10-40 mg, verapamil",
      "40-180 mg"
    ),
    regions = "literature data (healthy-volunteer DDI studies)",
    notes = paste(
      "Bottom-up prediction model, not a population fit. Transporter and",
      "metabolism parameters were taken from in vitro data and prior",
      "models (Table 1); physiology from Table S1; tissue:plasma ratios",
      "from Table S2. Clinical plasma data were used only for comparison.",
      "The paper's visual predictive check estimated variances on",
      "Vmax,OCT3, CLint,up,OCT1, CLint,up,OCT2, CLint,eff,MATE and the",
      "gastrointestinal transit rate but did not report them, so the model",
      "is typical-value only."
    )
  )

  # Per-perpetrator inputs for the parameterised perpetrator slot, one
  # column per drug. Parameter names match the ini() stems; to simulate a
  # different interaction, load a column into the model's p_* / ki_*_um /
  # kp_p / tissopt / cellopt parameters. Values are from Table 1 (fub, Rb,
  # CLint,met, CLint,OCT1/OCT2/MATE, CLrenal, Ki) and Table S2 (tissue:blood
  # ratios, divided by Rb). Molecular weights are from PubChem (non-paper
  # provenance). Three as-run choices follow the answered extraction
  # questions: p_clint_pd = 0 for the cimetidine transporter arm
  # (insensitive); rabeprazole p_clmet = 955.6 / fub = 23890 so that the
  # model's fub*CLint term applies CLint as a blood clearance without fub;
  # and the well-stirred p_clint_renal is back-calculated from CLrenal via
  # Eq 15. Ki,OCT3 > 100 uM (pyrimethamine) and the absent verapamil
  # Ki,MATE ship as 1e6 (no inhibition). cimetidine uses the transporter
  # structure (p_tr = 1); the other five use well-stirred (p_tr = 0).
  drugLibrary <- data.frame(
    parameter = c(
      "p_tr",
      "p_fub",
      "p_rbp",
      "p_clmet",
      "p_ka_duodenum",
      "p_ka_jejunum",
      "p_ka_ileum",
      "p_clint_up_oct1",
      "p_clint_up_oct2",
      "p_clint_eff_mate",
      "p_clint_renal",
      "p_clint_pd",
      "p_mw",
      "ki_oct1_um",
      "ki_oct2_um",
      "ki_oct3_um",
      "ki_mate_um",
      "tissopt",
      "cellopt",
      "kgb_p",
      "ksb_p",
      "kspb_p",
      "khb_p",
      "krb_p",
      "klub_p",
      "kheb_p",
      "kbrb_p",
      "kmub_p",
      "kadb_p",
      "kskb_p",
      "kreb_p"
    ),
    cimetidine = c(
      1,
      0.82,
      0.97,
      188.3,
      0.012,
      0.01472393,
      0.01655172,
      200,
      539.8,
      532.2,
      0,
      0,
      252.34,
      101,
      2.97,
      45.7,
      0.65,
      1,
      1,
      1.216495,
      1.216495,
      2.082474,
      1.505155,
      1.628866,
      2.42268,
      1.463918,
      1.670103,
      1.216495,
      1.948454,
      2.381443,
      0.01030928
    ),
    pyrimethamine = c(
      0,
      0.15,
      0.86,
      27.29,
      0.062,
      0.062,
      0.062,
      0,
      0,
      0,
      3.134521,
      0,
      248.71,
      4.46,
      0.61,
      1e6,
      0.02,
      1,
      1,
      3.604651,
      3.604651,
      4.186047,
      8.72093,
      7.906977,
      2.790698,
      3.139535,
      2.209302,
      1.430233,
      7.906977,
      2.44186,
      0.01162791
    ),
    trimethoprim = c(
      0,
      0.43,
      1.29,
      52.4,
      0.0355,
      0.0355,
      0.0355,
      0,
      0,
      0,
      193.2009,
      0,
      290.32,
      27.7,
      19.8,
      12.3,
      0.51,
      1,
      1,
      1.007752,
      1.007752,
      1.465116,
      1.472868,
      1.472868,
      1.581395,
      1.03876,
      1.186047,
      0.7596899,
      2.945736,
      1.503876,
      0.7751938
    ),
    ondansetron = c(
      0,
      0.325,
      0.83,
      1365.14,
      0.012,
      0.01472393,
      0.01655172,
      0,
      0,
      0,
      0,
      0,
      293.37,
      0.27,
      0.89,
      17.4,
      0.01,
      1,
      1,
      2.614458,
      0.8915663,
      0.8795181,
      1.409639,
      1.349398,
      1.783133,
      1.277108,
      2.39759,
      0.8915663,
      47.01205,
      3.759036,
      0.01204819
    ),
    rabeprazole = c(
      0,
      0.04,
      0.75,
      23890,
      0.029,
      0.029,
      0.029,
      0,
      0,
      0,
      0,
      0,
      359.44,
      3.0,
      5.7,
      3.0,
      4.60,
      1,
      1,
      0.2266667,
      0.1066667,
      0.1466667,
      0.1333333,
      0.1866667,
      0.2933333,
      0.2266667,
      0.08,
      0.1066667,
      0.1066667,
      0.3866667,
      0.01333333
    ),
    verapamil = c(
      0,
      0.14,
      0.68,
      17883,
      0.0156,
      0.0191411,
      0.02151724,
      0,
      0,
      0,
      0,
      0,
      454.60,
      9.62,
      3.24,
      3.6,
      1e6,
      1,
      1,
      32.16176,
      32.16176,
      36.85294,
      89.04412,
      80.44118,
      20.94118,
      28.33824,
      12.77941,
      11.14706,
      19.04412,
      18.36765,
      0.01470588
    ),
    stringsAsFactors = FALSE
  )

  ini({
    # ================= Physiology: 70-kg human (Table S1) ================
    # Gastrointestinal transit rate constants (1/min), Table S1.
    kt0 <- fixed(0.0462); label("Gastric emptying rate constant, stomach (1/min)") # Table S1 stomach
    kt1 <- fixed(0.0462); label("Duodenum transit rate constant (1/min)")          # Table S1 duodenum
    kt2 <- fixed(0.012);  label("Jejunum transit rate constant (1/min)")           # Table S1 jejunum
    kt3 <- fixed(0.0058); label("Ileum transit rate constant (1/min)")             # Table S1 ileum
    kt4 <- fixed(0.004);  label("Cecum transit rate constant (1/min)")             # Table S1 cecum
    kt5 <- fixed(0.0013); label("Colon transit rate constant (1/min)")             # Table S1 colon

    # Gut-wall (enterocyte) volumes (mL); gw_stomach = Table S1 stomach.
    v_gw_stomach  <- fixed(160);  label("Stomach wall volume (mL)")  # Table S1 stomach V
    v_gw_duodenum <- fixed(70);   label("Duodenum wall volume (mL)") # Table S1 duodenum V
    v_gw_jejunum  <- fixed(209);  label("Jejunum wall volume (mL)")  # Table S1 jejunum V
    v_gw_ileum    <- fixed(139);  label("Ileum wall volume (mL)")    # Table S1 ileum V
    v_gw_cecum    <- fixed(116);  label("Cecum wall volume (mL)")    # Table S1 cecum V
    v_gw_colon    <- fixed(1116); label("Colon wall volume (mL)")    # Table S1 colon V

    # Intestinal lumen volumes (mL), Eq 4 (duodenum, jejunum, ileum).
    vl_duodenum <- fixed(314.2);   label("Duodenum lumen volume (mL)") # Eq 4 Vi, duodenum
    vl_jejunum  <- fixed(2170.19); label("Jejunum lumen volume (mL)")  # Eq 4 Vi, jejunum
    vl_ileum    <- fixed(2609.05); label("Ileum lumen volume (mL)")    # Eq 4 Vi, ileum

    # Perfused-organ volumes (mL), Table S1.
    v_spleen <- fixed(190);   label("Spleen volume (mL)")  # Table S1 spleen
    v_hb     <- fixed(608.3); label("Hepatic blood volume (mL)")  # Eqs 7-8 Vh,b
    v_hc     <- fixed(1081.7); label("Hepatocyte volume (mL)")    # Eqs 7-8 Vh,c
    v_liver  <- fixed(1690);  label("Liver volume, well-stirred perpetrator (mL)") # Table S1 liver
    v_rb     <- fixed(33.8);  label("Renal blood volume (mL)")    # Eqs 12-13 Vr,b
    v_re     <- fixed(246.2); label("Renal tubule volume (mL)")   # Eqs 12-13 Vr,e
    v_kidney <- fixed(280);   label("Kidney volume, well-stirred perpetrator (mL)") # Table S1 kidney
    v_lung   <- fixed(1170);  label("Lung volume (mL)")     # Table S1 lung
    v_art    <- fixed(1730);  label("Arterial blood volume (mL)") # Table S1 artery
    v_ven    <- fixed(3470);  label("Venous blood volume (mL)")   # Table S1 venous
    v_heart  <- fixed(310);   label("Heart volume (mL)")    # Table S1 heart
    v_brain  <- fixed(1450);  label("Brain volume (mL)")    # Table S1 brain
    v_muscle <- fixed(35000); label("Muscle volume (mL)")   # Table S1 muscle
    v_adipose <- fixed(10000); label("Adipose volume (mL)") # Table S1 adipose
    v_skin   <- fixed(7800);  label("Skin volume (mL)")     # Table S1 skin
    v_other  <- fixed(5100);  label("Rest-of-body volume (mL)") # Table S1 rest of body

    # Blood flows (mL/min), Table S1.
    q_gw_stomach  <- fixed(38);  label("Stomach wall blood flow (mL/min)")  # Table S1 stomach Q
    q_gw_duodenum <- fixed(118); label("Duodenum wall blood flow (mL/min)") # Table S1 duodenum Q
    q_gw_jejunum  <- fixed(413); label("Jejunum wall blood flow (mL/min)")  # Table S1 jejunum Q
    q_gw_ileum    <- fixed(244); label("Ileum wall blood flow (mL/min)")    # Table S1 ileum Q
    q_gw_cecum    <- fixed(44);  label("Cecum wall blood flow (mL/min)")    # Table S1 cecum Q
    q_gw_colon    <- fixed(281); label("Colon wall blood flow (mL/min)")    # Table S1 colon Q
    q_spleen <- fixed(80);   label("Spleen blood flow (mL/min)")  # Table S1 spleen Q
    q_h      <- fixed(300);  label("Hepatic-artery blood flow (mL/min)") # Table S1 liver Q (hepatic artery)
    q_liver_out      <- fixed(1518); label("Total hepatic (outflow) blood flow (mL/min)") # Qh + Qsp + sum(Qgwi)
    q_r      <- fixed(1240); label("Renal blood flow (mL/min)")   # Table S1 kidney Q
    q_tot    <- fixed(5600); label("Cardiac output (mL/min)")     # Table S1 lung Q
    q_heart  <- fixed(240);  label("Heart blood flow (mL/min)")   # Table S1 heart Q
    q_brain  <- fixed(700);  label("Brain blood flow (mL/min)")   # Table S1 brain Q
    q_muscle <- fixed(750);  label("Muscle blood flow (mL/min)")  # Table S1 muscle Q
    q_adipose <- fixed(260); label("Adipose blood flow (mL/min)") # Table S1 adipose Q
    q_skin   <- fixed(300);  label("Skin blood flow (mL/min)")    # Table S1 skin Q
    q_other  <- fixed(592);  label("Rest-of-body blood flow (mL/min)") # Table S1 rest of body Q

    # ================= Metformin absorption (Eqs 2-3) ===================
    peff_met   <- fixed(0.00311); label("Metformin effective permeability Peff (cm/min)") # Table 1 Met Peff
    r_duodenum <- fixed(2.0);     label("Duodenum radius (cm)")  # Eq 3 ri, duodenum
    r_jejunum  <- fixed(1.63);    label("Jejunum radius (cm)")   # Eq 3 ri, jejunum
    r_ileum    <- fixed(1.45);    label("Ileum radius (cm)")     # Eq 3 ri, ileum
    pbsf_duodenum <- fixed(2790.65);  label("Duodenum S9 protein (mg)") # Eq 2 PBSF, duodenum
    pbsf_jejunum  <- fixed(14465.81); label("Jejunum S9 protein (mg)")  # Eq 2 PBSF, jejunum
    pbsf_ileum    <- fixed(11219.5);  label("Ileum S9 protein (mg)")    # Eq 2 PBSF, ileum
    # OCT3 / PMAT Michaelis-Menten constants, Table 1 / Methods. Shipped
    # AS-RUN (q2): Vmax,OCT3 12.08 and Vmax,PMAT 15.28 nmol/mg/min are used
    # as 12.08/1000 and 15.28/1000, and Km,OCT3 2.46 mM and Km,PMAT 1.68 mM
    # are used as the bare numbers against the ug/mL lumen concentration,
    # with no molar<->mass conversion. This dimensionally inconsistent form
    # is the one that reproduces Table S4 and the knock-outs; see Errata.
    vmax_oct3 <- fixed(0.01208); label("OCT3 maximum velocity, as-run (Vmax 12.08 nmol/mg/min / 1000)") # Methods 'Vmax,OCT3 12.08 nmol/mg/min'
    km_oct3   <- fixed(2.46);    label("OCT3 Michaelis constant, as-run (Km 2.46 mM used unconverted)")  # Methods 'Km,OCT3 2.46 mM'
    vmax_pmat <- fixed(0.01528); label("PMAT maximum velocity, as-run (Vmax 15.28 nmol/mg/min / 1000)")  # Methods 'Vmax,PMAT 15.28 nmol/mg/min'
    km_pmat   <- fixed(1.68);    label("PMAT Michaelis constant, as-run (Km 1.68 mM used unconverted)")   # Methods 'Km,PMAT 1.68 mM'
    # Regional OCT3 expression scaling Tsf,ex,i (jejunum = 1).
    tex_duodenum <- fixed(1.19);  label("OCT3 expression scaling, duodenum (unitless)") # Methods 'Tsf,ex 1.19:1:1.155'
    tex_jejunum  <- fixed(1.0);   label("OCT3 expression scaling, jejunum (unitless)")  # Methods Tsf,ex jejunum = 1
    tex_ileum    <- fixed(1.155); label("OCT3 expression scaling, ileum (unitless)")    # Methods 'Tsf,ex ... 1.155'
    # pH-dependent activity scaling Tsf,pH,i for OCT3 and PMAT.
    tph_oct3_duodenum <- fixed(0.35); label("OCT3 pH scaling, duodenum (unitless)") # Methods 'Tsf,pH OCT3 0.35, 0.74, 1'
    tph_oct3_jejunum  <- fixed(0.74); label("OCT3 pH scaling, jejunum (unitless)")  # Methods Tsf,pH OCT3 jejunum
    tph_oct3_ileum    <- fixed(1.0);  label("OCT3 pH scaling, ileum (unitless)")    # Methods Tsf,pH OCT3 ileum
    tph_pmat_duodenum <- fixed(1.0);  label("PMAT pH scaling, duodenum (unitless)") # Methods 'Tsf,pH PMAT 1, 0.48, 0.11'
    tph_pmat_jejunum  <- fixed(0.48); label("PMAT pH scaling, jejunum (unitless)")  # Methods Tsf,pH PMAT jejunum
    tph_pmat_ileum    <- fixed(0.11); label("PMAT pH scaling, ileum (unitless)")    # Methods Tsf,pH PMAT ileum

    # ================= Metformin disposition (Table 1 Met) ==============
    fub_met <- fixed(1);    label("Metformin unbound fraction in blood (fraction)") # Table 1 Met fub
    rbp_met <- fixed(1);    label("Metformin blood-to-plasma ratio Rb (unitless)")  # Table 1 Met Rb
    clint_up_oct1 <- fixed(63.95);  label("Hepatic OCT1 uptake intrinsic clearance (mL/min)")  # Table 1 Met CLint,OCT1
    clint_met     <- fixed(110.57); label("Hepatic metabolic clearance (mL/min)")              # Table 1 Met CLint,met,h
    clint_up_oct2 <- fixed(256.45); label("Renal OCT2 uptake intrinsic clearance (mL/min)")    # Table 1 Met CLint,OCT2
    clint_eff_mate <- fixed(299.8); label("Renal MATE efflux intrinsic clearance (mL/min)")    # Table 1 Met CLint,MATE
    raf  <- fixed(3.0);   label("Renal relative activity factor RAF (unitless)") # Methods 'RAF ... 3.0'
    gfr  <- fixed(132.9); label("Glomerular filtration rate (mL/min)")           # Table 1 Met CLrenal / GFR 132.9
    # Hepatocyte<->blood passive-diffusion clearance CLint,pd (Eqs 7-8),
    # printed nowhere. Back-solved (q1): six independent printed outputs
    # (control AUC, no-OCT1 +15.67%, no-OCT2 +160.30%, no-OCT2 AUC 12.50,
    # AUC500 7.68, AUC1000 12.55) all cross between 10 and 12 mL/min.
    clint_pd <- fixed(11); label("Hepatic passive-diffusion clearance, back-solved (mL/min)") # back-solved from Fig 6 / sensitivity-text AUCs and Table S4

    # Metformin tissue:blood ratios = Table S2 Met column / Rb (Rb = 1).
    kgb_met  <- fixed(1.5);  label("Metformin gut-wall:blood ratio (unitless)")  # Table S2 Met intestine 1.5 / Rb 1
    ksb_met  <- fixed(1.5);  label("Metformin stomach-wall:blood ratio (unitless)") # Table S2 Met stomach 1.5
    kspb_met <- fixed(1.66); label("Metformin spleen:blood ratio (unitless)")    # Table S2 Met spleen 1.66
    khb_met  <- fixed(1.6);  label("Metformin liver:blood ratio (unitless)")     # Table S2 Met liver 1.6
    krb_met  <- fixed(1.49); label("Metformin kidney:blood ratio (unitless)")    # Table S2 Met kidney 1.49
    klub_met <- fixed(1.46); label("Metformin lung:blood ratio (unitless)")      # Table S2 Met lung 1.46
    kheb_met <- fixed(1.47); label("Metformin heart:blood ratio (unitless)")     # Table S2 Met heart 1.47
    kbrb_met <- fixed(1.72); label("Metformin brain:blood ratio (unitless)")     # Table S2 Met brain 1.72
    kmub_met <- fixed(1.7);  label("Metformin muscle:blood ratio (unitless)")    # Table S2 Met muscle 1.7
    kadb_met <- fixed(0.18); label("Metformin adipose:blood ratio (unitless)")   # Table S2 Met adipose 0.18
    kskb_met <- fixed(1.11); label("Metformin skin:blood ratio (unitless)")      # Table S2 Met skin 1.11
    kreb_met <- fixed(1.2);  label("Metformin rest-of-body:blood ratio (unitless)") # Table S2 Met rest 1.2

    # ================= Perpetrator slot (ships = cimetidine) ============
    # tr = 1 selects the transporter-resolved liver/kidney structure
    # (cimetidine); tr = 0 selects well-stirred liver and kidney (the five
    # other perpetrators). Load another drug from the drugLibrary metadata.
    p_tr <- fixed(1); label("Perpetrator structure switch: 1 = transporter, 0 = well-stirred") # cimetidine uses transporter
    p_fub <- fixed(0.82); label("Perpetrator unbound fraction in blood (fraction)") # Table 1 Cim fub
    p_rbp <- fixed(0.97); label("Perpetrator blood-to-plasma ratio Rb (unitless)")  # Table 1 Cim Rb
    p_clmet <- fixed(188.3); label("Perpetrator hepatic metabolic clearance (mL/min)") # Table 1 Cim CLint,met,h
    # Perpetrator segmental absorption rate constants (1/min). Cimetidine:
    # ka,i = 2 * Peff / r_i with Peff 0.012 cm/min (Table 1 Cim Peff).
    p_ka_duodenum <- fixed(0.012);      label("Perpetrator duodenum absorption rate constant (1/min)") # 2*0.012/2.0 (Cim Peff)
    p_ka_jejunum  <- fixed(0.01472393); label("Perpetrator jejunum absorption rate constant (1/min)")  # 2*0.012/1.63
    p_ka_ileum    <- fixed(0.01655172); label("Perpetrator ileum absorption rate constant (1/min)")    # 2*0.012/1.45
    # Transporter-arm clearances (cimetidine); zero for well-stirred drugs.
    p_clint_up_oct1 <- fixed(200);   label("Perpetrator hepatic OCT1 uptake clearance (mL/min)") # Table 1 Cim CLint,OCT1
    p_clint_up_oct2 <- fixed(539.8); label("Perpetrator renal OCT2 uptake clearance (mL/min)")   # Table 1 Cim CLint,OCT2
    p_clint_eff_mate <- fixed(532.2); label("Perpetrator renal MATE efflux clearance (mL/min)")  # Table 1 Cim CLint,MATE
    # Well-stirred renal intrinsic clearance (zero for cimetidine);
    # back-calculated from CLrenal via Eq 15 for the five other drugs.
    p_clint_renal <- fixed(0); label("Perpetrator well-stirred renal intrinsic clearance (mL/min)") # cimetidine: renal via transporters
    # Hepatic passive-diffusion clearance for the cimetidine transporter
    # arm; cimetidine is insensitive to it (q1), so 0.
    p_clint_pd <- fixed(0); label("Perpetrator hepatic passive-diffusion clearance (mL/min)") # cimetidine insensitive (q1)
    p_mw <- fixed(252.34); label("Perpetrator molecular weight (g/mol)") # cimetidine MW (PubChem CID 2756)
    # Perpetrator inhibition constants (uM), Table 1. Converted to ug/mL in
    # model() by Ki[uM] * MW / 1000. >100 uM OCT3 (pyrimethamine) ships as
    # 1e6 (no inhibition).
    ki_oct1_um <- fixed(101);  label("Perpetrator OCT1 inhibition constant (uM)")  # Table 1 Cim Ki,OCT1
    ki_oct2_um <- fixed(2.97); label("Perpetrator OCT2 inhibition constant (uM)")  # Table 1 Cim Ki,OCT2
    ki_oct3_um <- fixed(45.7); label("Perpetrator OCT3 inhibition constant (uM)")  # Table 1 Cim Ki,OCT3
    ki_mate_um <- fixed(0.65); label("Perpetrator MATE inhibition constant (uM)")  # Table 1 Cim Ki,MATE
    # Inhibitor-concentration selectors for the well-stirred perpetrator
    # arm (unused when tr = 1): 1 = use the well-stirred tissue
    # concentration A/V; 0 = use the blood-equilibrium A/V/Kp.
    tissopt <- fixed(1); label("Well-stirred liver inhibitor-concentration selector (1 = tissue)") # fub*A/V form
    cellopt <- fixed(1); label("Well-stirred tubule inhibitor-concentration selector (1 = tissue)") # fub*A/V form

    # Perpetrator tissue:blood ratios = Table S2 Cim column / Rb.
    kgb_p  <- fixed(1.216495); label("Perpetrator gut-wall:blood ratio (unitless)")  # Table S2 Cim intestine 1.18 / Rb 0.97
    ksb_p  <- fixed(1.216495); label("Perpetrator stomach-wall:blood ratio (unitless)") # Table S2 Cim stomach 1.18 / 0.97
    kspb_p <- fixed(2.082474); label("Perpetrator spleen:blood ratio (unitless)")    # Table S2 Cim spleen 2.02 / 0.97
    khb_p  <- fixed(1.505155); label("Perpetrator liver:blood ratio (unitless)")     # Table S2 Cim liver 1.46 / 0.97
    krb_p  <- fixed(1.628866); label("Perpetrator kidney:blood ratio (unitless)")    # Table S2 Cim kidney 1.58 / 0.97
    klub_p <- fixed(2.42268);  label("Perpetrator lung:blood ratio (unitless)")      # Table S2 Cim lung 2.35 / 0.97
    kheb_p <- fixed(1.463918); label("Perpetrator heart:blood ratio (unitless)")     # Table S2 Cim heart 1.42 / 0.97
    kbrb_p <- fixed(1.670103); label("Perpetrator brain:blood ratio (unitless)")     # Table S2 Cim brain 1.62 / 0.97
    kmub_p <- fixed(1.216495); label("Perpetrator muscle:blood ratio (unitless)")    # Table S2 Cim muscle 1.18 / 0.97
    kadb_p <- fixed(1.948454); label("Perpetrator adipose:blood ratio (unitless)")   # Table S2 Cim adipose 1.89 / 0.97
    kskb_p <- fixed(2.381443); label("Perpetrator skin:blood ratio (unitless)")      # Table S2 Cim skin 2.31 / 0.97
    kreb_p <- fixed(0.01030928); label("Perpetrator rest-of-body:blood ratio (unitless)") # Table S2 Cim rest 0.01 / 0.97

    # ================= Residual error ===================================
    # Deterministic prediction model; the paper reports no residual-error
    # estimate, so it is fixed at zero rather than invented.
    propSd <- fixed(0); label("Proportional residual error (fraction; ZERO - not reported)") # not reported (VPC variances unpublished)
  })

  model({
    # ================= Derived rate constants ===========================
    # Metformin segmental absorption, Eq 3 (ka,i = 2 Peff / r_i).
    ka_duodenum <- 2 * peff_met / r_duodenum
    ka_jejunum  <- 2 * peff_met / r_jejunum
    ka_ileum    <- 2 * peff_met / r_ileum

    # Perpetrator inhibition constants in ug/mL (Ki[uM] * MW / 1000).
    ki_oct1 <- ki_oct1_um * p_mw / 1000
    ki_oct2 <- ki_oct2_um * p_mw / 1000
    ki_oct3 <- ki_oct3_um * p_mw / 1000
    ki_mate <- ki_mate_um * p_mw / 1000

    # ================= Perpetrator gut lumen (Eq 5) =====================
    d/dt(stomach_perpetrator)  <- -kt0 * stomach_perpetrator
    pj_duo <- p_ka_duodenum * duodenum_perpetrator
    pj_jej <- p_ka_jejunum * jejunum_perpetrator
    pj_ile <- p_ka_ileum * ileum_perpetrator
    d/dt(duodenum_perpetrator) <- kt0 * stomach_perpetrator - kt1 * duodenum_perpetrator - pj_duo
    d/dt(jejunum_perpetrator)  <- kt1 * duodenum_perpetrator - kt2 * jejunum_perpetrator - pj_jej
    d/dt(ileum_perpetrator)    <- kt2 * jejunum_perpetrator - kt3 * ileum_perpetrator - pj_ile
    d/dt(cecum_perpetrator)    <- kt3 * ileum_perpetrator - kt4 * cecum_perpetrator
    d/dt(colon_perpetrator)    <- kt4 * cecum_perpetrator - kt5 * colon_perpetrator

    # ================= Perpetrator gut wall (Eq 6) ======================
    p_cart <- arterial_perpetrator / v_art
    d/dt(gw_stomach_perpetrator)  <- q_gw_stomach * (p_cart - gw_stomach_perpetrator / v_gw_stomach / ksb_p)
    d/dt(gw_duodenum_perpetrator) <- q_gw_duodenum * (p_cart - gw_duodenum_perpetrator / v_gw_duodenum / kgb_p) + pj_duo
    d/dt(gw_jejunum_perpetrator)  <- q_gw_jejunum * (p_cart - gw_jejunum_perpetrator / v_gw_jejunum / kgb_p) + pj_jej
    d/dt(gw_ileum_perpetrator)    <- q_gw_ileum * (p_cart - gw_ileum_perpetrator / v_gw_ileum / kgb_p) + pj_ile
    d/dt(gw_cecum_perpetrator)    <- q_gw_cecum * (p_cart - gw_cecum_perpetrator / v_gw_cecum / kgb_p)
    d/dt(gw_colon_perpetrator)    <- q_gw_colon * (p_cart - gw_colon_perpetrator / v_gw_colon / kgb_p)
    d/dt(spleen_perpetrator)      <- q_spleen * (p_cart - spleen_perpetrator / v_spleen / kspb_p)
    p_portal <- q_gw_stomach * gw_stomach_perpetrator / v_gw_stomach / ksb_p +
      q_gw_duodenum * gw_duodenum_perpetrator / v_gw_duodenum / kgb_p +
      q_gw_jejunum * gw_jejunum_perpetrator / v_gw_jejunum / kgb_p +
      q_gw_ileum * gw_ileum_perpetrator / v_gw_ileum / kgb_p +
      q_gw_cecum * gw_cecum_perpetrator / v_gw_cecum / kgb_p +
      q_gw_colon * gw_colon_perpetrator / v_gw_colon / kgb_p +
      q_spleen * spleen_perpetrator / v_spleen / kspb_p

    # ================= Perpetrator liver: transporter arm (tr = 1) ======
    p_chb <- liver_blood_perpetrator / v_hb
    p_chc <- hepatocyte_perpetrator / v_hc / khb_p
    p_up1 <- p_clint_up_oct1 * p_fub * p_chb
    p_pd  <- p_fub * p_clint_pd * (p_chc - p_chb)
    d/dt(liver_blood_perpetrator) <- p_tr * (q_h * p_cart + p_portal - q_liver_out * p_chb - p_up1 + p_pd)
    d/dt(hepatocyte_perpetrator)  <- p_tr * (p_up1 - p_pd - p_fub * p_clmet * p_chc)

    # ================= Perpetrator liver: well-stirred arm (tr = 0) ======
    p_ch <- liver_perpetrator / v_liver / khb_p
    d/dt(liver_perpetrator) <- (1 - p_tr) * (q_h * p_cart + p_portal - q_liver_out * p_ch - p_fub * p_clmet * p_ch)

    # ================= Perpetrator kidney: transporter arm (tr = 1) =====
    p_crb <- kidney_blood_perpetrator / v_rb
    p_cre <- tubule_perpetrator / v_re / krb_p
    d/dt(kidney_blood_perpetrator) <- p_tr * (q_r * (p_cart - p_crb) -
      raf * p_clint_up_oct2 * p_fub * p_crb + p_clint_up_oct2 * p_fub * p_cre - p_fub * gfr * p_crb)
    d/dt(tubule_perpetrator) <- p_tr * (raf * p_clint_up_oct2 * p_fub * p_crb -
      p_clint_up_oct2 * p_fub * p_cre - raf * p_clint_eff_mate * p_fub * p_cre)

    # ================= Perpetrator kidney: well-stirred arm (tr = 0) =====
    p_ck <- kidney_perpetrator / v_kidney / krb_p
    d/dt(kidney_perpetrator) <- (1 - p_tr) * (q_r * (p_cart - p_ck) - p_fub * p_clint_renal * p_ck)

    # ================= Perpetrator blood and tissues ====================
    p_cven <- venous_perpetrator / v_ven
    d/dt(lung_perpetrator) <- q_tot * (p_cven - lung_perpetrator / v_lung / klub_p)
    d/dt(arterial_perpetrator) <- q_tot * (lung_perpetrator / v_lung / klub_p - p_cart)
    d/dt(heart_perpetrator)   <- q_heart * (p_cart - heart_perpetrator / v_heart / kheb_p)
    d/dt(brain_perpetrator)   <- q_brain * (p_cart - brain_perpetrator / v_brain / kbrb_p)
    d/dt(muscle_perpetrator)  <- q_muscle * (p_cart - muscle_perpetrator / v_muscle / kmub_p)
    d/dt(adipose_perpetrator) <- q_adipose * (p_cart - adipose_perpetrator / v_adipose / kadb_p)
    d/dt(skin_perpetrator)    <- q_skin * (p_cart - skin_perpetrator / v_skin / kskb_p)
    d/dt(other_perpetrator)   <- q_other * (p_cart - other_perpetrator / v_other / kreb_p)
    d/dt(venous_perpetrator)  <- q_heart * heart_perpetrator / v_heart / kheb_p +
      q_brain * brain_perpetrator / v_brain / kbrb_p +
      q_muscle * muscle_perpetrator / v_muscle / kmub_p +
      q_adipose * adipose_perpetrator / v_adipose / kadb_p +
      q_skin * skin_perpetrator / v_skin / kskb_p +
      q_other * other_perpetrator / v_other / kreb_p +
      p_tr * q_liver_out * p_chb + (1 - p_tr) * q_liver_out * p_ch +
      p_tr * q_r * p_crb + (1 - p_tr) * q_r * p_ck - q_tot * p_cven

    # ================= Perpetrator inhibitor concentrations =============
    # The inhibitor concentration presented to each metformin transporter.
    # Transporter arm: the resolved blood (liver_blood, kidney_blood) and
    # tubule amounts. Well-stirred arm: the well-stirred tissue amount,
    # selected by tissopt / cellopt (tissue A/V vs blood-equilibrium
    # A/V/Kp). Multiplied by p_fub in the metformin equations.
    ci_h  <- p_tr * liver_blood_perpetrator / v_hb +
      (1 - p_tr) * (tissopt * liver_perpetrator / v_liver + (1 - tissopt) * liver_perpetrator / v_liver / khb_p)
    ci_rb <- p_tr * kidney_blood_perpetrator / v_rb +
      (1 - p_tr) * (tissopt * kidney_perpetrator / v_kidney + (1 - tissopt) * kidney_perpetrator / v_kidney / krb_p)
    ci_re <- p_tr * tubule_perpetrator / v_re +
      (1 - p_tr) * (cellopt * kidney_perpetrator / v_kidney + (1 - cellopt) * kidney_perpetrator / v_kidney / krb_p)

    # ================= Metformin gut lumen (Eqs 1-2) ====================
    # M-M driven by the lumen concentration A_i / V_i (q2), not the printed
    # enterocyte concentration. Intestinal OCT3 is competitively inhibited
    # by the perpetrator lumen concentration (Eq 4).
    d/dt(stomach) <- -kt0 * stomach
    c_duo <- duodenum / vl_duodenum
    c_jej <- jejunum / vl_jejunum
    c_ile <- ileum / vl_ileum
    j_duo <- 0.3 * ka_duodenum * duodenum +
      pbsf_duodenum * vmax_oct3 / (1 + duodenum_perpetrator / vl_duodenum / ki_oct3) * c_duo / (km_oct3 + c_duo) * tex_duodenum * tph_oct3_duodenum +
      pbsf_duodenum * vmax_pmat * c_duo / (km_pmat + c_duo) * tph_pmat_duodenum
    j_jej <- 0.3 * ka_jejunum * jejunum +
      pbsf_jejunum * vmax_oct3 / (1 + jejunum_perpetrator / vl_jejunum / ki_oct3) * c_jej / (km_oct3 + c_jej) * tex_jejunum * tph_oct3_jejunum +
      pbsf_jejunum * vmax_pmat * c_jej / (km_pmat + c_jej) * tph_pmat_jejunum
    j_ile <- 0.3 * ka_ileum * ileum +
      pbsf_ileum * vmax_oct3 / (1 + ileum_perpetrator / vl_ileum / ki_oct3) * c_ile / (km_oct3 + c_ile) * tex_ileum * tph_oct3_ileum +
      pbsf_ileum * vmax_pmat * c_ile / (km_pmat + c_ile) * tph_pmat_ileum
    d/dt(duodenum) <- kt0 * stomach - kt1 * duodenum - j_duo
    d/dt(jejunum)  <- kt1 * duodenum - kt2 * jejunum - j_jej
    d/dt(ileum)    <- kt2 * jejunum - kt3 * ileum - j_ile
    d/dt(cecum)    <- kt3 * ileum - kt4 * cecum
    d/dt(colon)    <- kt4 * cecum - kt5 * colon
    d/dt(a_feces)  <- kt5 * colon

    # ================= Metformin gut wall (Eq 3) ========================
    cart <- arterial / v_art
    d/dt(gw_stomach)  <- q_gw_stomach * (cart - gw_stomach / v_gw_stomach / ksb_met)
    d/dt(gw_duodenum) <- q_gw_duodenum * (cart - gw_duodenum / v_gw_duodenum / kgb_met) + j_duo
    d/dt(gw_jejunum)  <- q_gw_jejunum * (cart - gw_jejunum / v_gw_jejunum / kgb_met) + j_jej
    d/dt(gw_ileum)    <- q_gw_ileum * (cart - gw_ileum / v_gw_ileum / kgb_met) + j_ile
    d/dt(gw_cecum)    <- q_gw_cecum * (cart - gw_cecum / v_gw_cecum / kgb_met)
    d/dt(gw_colon)    <- q_gw_colon * (cart - gw_colon / v_gw_colon / kgb_met)
    d/dt(spleen)      <- q_spleen * (cart - spleen / v_spleen / kspb_met)

    # ================= Metformin liver (Eqs 7-8) ========================
    # OCT1 uptake competitively inhibited by the perpetrator (p_fub*ci_h).
    chb <- liver_blood / v_hb
    chc <- hepatocyte / v_hc / khb_met
    up1 <- clint_up_oct1 / (1 + p_fub * ci_h / ki_oct1) * fub_met * chb
    pd  <- fub_met * clint_pd * (chc - chb)
    d/dt(liver_blood) <- q_h * cart +
      q_gw_stomach * gw_stomach / v_gw_stomach / ksb_met +
      q_gw_duodenum * gw_duodenum / v_gw_duodenum / kgb_met +
      q_gw_jejunum * gw_jejunum / v_gw_jejunum / kgb_met +
      q_gw_ileum * gw_ileum / v_gw_ileum / kgb_met +
      q_gw_cecum * gw_cecum / v_gw_cecum / kgb_met +
      q_gw_colon * gw_colon / v_gw_colon / kgb_met +
      q_spleen * spleen / v_spleen / kspb_met - q_liver_out * chb - up1 + pd
    d/dt(hepatocyte) <- up1 - pd - fub_met * clint_met * chc

    # ================= Metformin kidney (Eqs 12-13) =====================
    # OCT2 uptake and back-flux, and MATE efflux, each competitively
    # inhibited by the perpetrator.
    crb <- kidney_blood / v_rb
    cre <- tubule / v_re / krb_met
    d/dt(kidney_blood) <- q_r * (cart - crb) -
      raf * clint_up_oct2 / (1 + p_fub * ci_rb / ki_oct2) * fub_met * crb +
      clint_up_oct2 / (1 + p_fub * ci_re / ki_oct2) * fub_met * cre -
      fub_met * gfr * crb
    d/dt(tubule) <- raf * clint_up_oct2 / (1 + p_fub * ci_rb / ki_oct2) * fub_met * crb -
      clint_up_oct2 / (1 + p_fub * ci_re / ki_oct2) * fub_met * cre -
      raf * clint_eff_mate / (1 + p_fub * ci_re / ki_mate) * fub_met * cre
    d/dt(a_urine) <- fub_met * gfr * crb +
      raf * clint_eff_mate / (1 + p_fub * ci_re / ki_mate) * fub_met * cre

    # ================= Metformin blood and tissues (Eqs 16-17) ==========
    cven <- venous / v_ven
    d/dt(lung) <- q_tot * (cven - lung / v_lung / klub_met)
    d/dt(arterial) <- q_tot * (lung / v_lung / klub_met - cart)
    d/dt(heart)   <- q_heart * (cart - heart / v_heart / kheb_met)
    d/dt(brain)   <- q_brain * (cart - brain / v_brain / kbrb_met)
    d/dt(muscle)  <- q_muscle * (cart - muscle / v_muscle / kmub_met)
    d/dt(adipose) <- q_adipose * (cart - adipose / v_adipose / kadb_met)
    d/dt(skin)    <- q_skin * (cart - skin / v_skin / kskb_met)
    d/dt(other)   <- q_other * (cart - other / v_other / kreb_met)
    d/dt(venous)  <- q_heart * heart / v_heart / kheb_met +
      q_brain * brain / v_brain / kbrb_met +
      q_muscle * muscle / v_muscle / kmub_met +
      q_adipose * adipose / v_adipose / kadb_met +
      q_skin * skin / v_skin / kskb_met +
      q_other * other / v_other / kreb_met +
      q_liver_out * chb + q_r * crb - q_tot * cven

    # ================= Outputs ==========================================
    Cc <- cven / rbp_met           # metformin plasma concentration (ug/mL)
    Cliver <- chc                  # metformin hepatocyte concentration (ug/mL)
    Cperp <- p_cven / p_rbp        # perpetrator plasma concentration (ug/mL)

    Cc ~ prop(propSd)
  })
}
