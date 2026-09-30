Xu_2020_clopidogrel_pbpk <- function() {
  description <- paste(
    "PBPK-PD (whole-body, perfusion-limited, Phoenix WinNonlin 8.2).",
    "Clopidogrel, 2-oxo-clopidogrel and the active thiol metabolite",
    "(CLOP-AM, H4) with inhibition of platelet aggregation (IPA) after oral",
    "clopidogrel in a 70-kg adult, with CYP2C19 phenotype (UM/EM/IM/PM),",
    "coronary artery disease (CAD) and diabetes mellitus (DM) covariates.",
    "Fourteen perfusion-limited tissue/blood compartments per analyte",
    "(lung, heart, spleen, liver, kidney, brain, adipose, muscle, skin, rest",
    "of body, stomach wall, venous and arterial blood) plus five gut-wall",
    "segments. Oral clopidogrel empties from the stomach lumen through a",
    "five-segment gut lumen; the duodenum, jejunum and ileum absorb it with",
    "P-gp efflux back to the lumen. Hepatic CYP1A2/2B6/2C19 Michaelis-Menten",
    "oxidation forms 2-oxo-clopidogrel, CYP2B6/2C9/3A4/2C19 oxidation of",
    "2-oxo-clopidogrel forms CLOP-AM, and CES1 hydrolysis inactivates all",
    "three analytes. Unbound venous CLOP-AM irreversibly inactivates",
    "platelets in a turnover model of normalized maximal platelet",
    "aggregation M, and IPA = (1 - M) * 100. CAD lowers all blood flows to",
    "0.9-fold and platelet responsiveness (kirre) to 0.7-fold; DM rescales",
    "CYP and CES1 activities and replaces the gastrointestinal transit",
    "rates. Deterministic typical-value model (the variances of the paper's",
    "visual predictive check were not reported)."
  )
  reference <- paste(
    "Xu RJ, Kong WM, An XF, Zou JJ, Liu L, Liu XD (2020).",
    "Physiologically-Based Pharmacokinetic-Pharmacodynamics Model",
    "Characterizing CYP2C19 Polymorphisms to Predict Clopidogrel",
    "Pharmacokinetics and Its Anti-Platelet Aggregation Effect Following",
    "Oral Administration to Coronary Artery Disease Patients With or",
    "Without Diabetes. Front Pharmacol 11:593982.",
    "doi:10.3389/fphar.2020.593982.",
    "Equations 1-18; physiology from Table 1; physicochemical parameters",
    "and Kt/p from Table 2; CYP kinetics from Table 3; all other values",
    "from the Methods text."
  )
  vignette <- "Xu_2020_clopidogrel_pbpk"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  paper_specific_compartments <- c(
    "wall_stomach",
    "wall_duodenum",
    "wall_jejunum",
    "wall_ileum",
    "wall_cecum",
    "wall_colon",
    "wall_stomach_oxoclop",
    "wall_duodenum_oxoclop",
    "wall_jejunum_oxoclop",
    "wall_ileum_oxoclop",
    "wall_cecum_oxoclop",
    "wall_colon_oxoclop",
    "wall_stomach_h4",
    "wall_duodenum_h4",
    "wall_jejunum_h4",
    "wall_ileum_h4",
    "wall_cecum_h4",
    "wall_colon_h4",
    "aggregation"
  )

  covariateData <- list(
    CYP2C19_UM = list(
      description = "CYP2C19 ultrarapid-metabolizer phenotype indicator (1 = UM: *1/*17 or *17/*17).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (extensive metabolizer *1/*1 when CYP2C19_IM and CYP2C19_PM are also 0)",
      notes = paste(
        "Multiplies the CYP2C19 Vmax of both oxidation steps by 1.58",
        "(Table 3 footnote b; 7.52 -> 11.88 and 9.06 -> 14.31",
        "pmol/pmol P450/min). CYP2C19_UM, CYP2C19_IM and CYP2C19_PM are",
        "mutually exclusive; all three 0 selects the EM phenotype."
      ),
      source_name = "UM"
    ),
    CYP2C19_IM = list(
      description = "CYP2C19 intermediate-metabolizer phenotype indicator (1 = IM: *1/*2, *2/*17 or *1/*3).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (extensive metabolizer *1/*1 when CYP2C19_UM and CYP2C19_PM are also 0)",
      notes = paste(
        "Multiplies the CYP2C19 Vmax of both oxidation steps by 0.5",
        "(Table 3 footnote c). Unlike Zhao 2018, the reference here is EM",
        "only; UM is flagged separately by CYP2C19_UM."
      ),
      source_name = "IM"
    ),
    CYP2C19_PM = list(
      description = "CYP2C19 poor-metabolizer phenotype indicator (1 = PM: *2/*2 or *2/*3).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (extensive metabolizer *1/*1 when CYP2C19_UM and CYP2C19_IM are also 0)",
      notes = "Sets the CYP2C19 Vmax of both oxidation steps to 0 (Table 3 footnote c).",
      source_name = "PM"
    ),
    DIS_IHD = list(
      description = "Coronary artery disease indicator (1 = CAD patient on aspirin, 0 = healthy).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (healthy volunteer)",
      notes = paste(
        "CAD scales every blood flow, including cardiac output, by",
        "Qtotal,CAD / Qtotal,health = 0.90 (Eq 18; Table 1 'CAD' column)",
        "and scales kirre by 0.7 (Methods 'PBPK-PD Model in CAD Without DM",
        "Patients'). The paper's DM population is CAD + DM, so a DM patient",
        "is DIS_IHD = 1 and DIS_DIAB = 1; DIS_DIAB = 1 with DIS_IHD = 0 is",
        "outside the paper's scenarios."
      ),
      source_name = "CAD"
    ),
    DIS_DIAB = list(
      description = "Diabetes mellitus indicator (1 = DM, 0 = no DM).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no diabetes)",
      notes = paste(
        "DM multiplies the CYP1A2, CYP2B6, CYP2C9, CYP2C19 and CYP3A4 Vmax",
        "by 1.23, 0.55, 1.26, 0.54 and 0.62 and the three CES1 intrinsic",
        "clearances by 1.27, and replaces the stomach..colon transit rate",
        "constants with 2.31, 2.30, 0.99, 1.32, 0.20 and 0.04 1/h (Methods",
        "'PBPK-PD Model in CAD Patients With DM'). kirre in DM equals the",
        "CAD value, so DM adds no further kirre change."
      ),
      source_name = "DM"
    )
  )

  # Every state is an amount in mg. The metabolite states hold the mass
  # formed 1:1 from the metabolized parent mass (Eqs 13 and 15 as run; no
  # molecular-weight ratio), so a 2-oxo-clopidogrel or CLOP-AM "mg" is the
  # paper's mass unit, not a stoichiometric one.
  compartmentData <- list(
    stomach = list(analyte = "clopidogrel", units = "mg", specimen = "administration site", verified = TRUE),
    duodenum = list(analyte = "clopidogrel", units = "mg", specimen = "administration site", verified = TRUE),
    jejunum = list(analyte = "clopidogrel", units = "mg", specimen = "administration site", verified = TRUE),
    ileum = list(analyte = "clopidogrel", units = "mg", specimen = "administration site", verified = TRUE),
    cecum = list(analyte = "clopidogrel", units = "mg", specimen = "administration site", verified = TRUE),
    colon = list(analyte = "clopidogrel", units = "mg", specimen = "administration site", verified = TRUE),
    wall_stomach = list(analyte = "clopidogrel", units = "mg", specimen = "tissue", verified = TRUE),
    wall_duodenum = list(analyte = "clopidogrel", units = "mg", specimen = "tissue", verified = TRUE),
    wall_jejunum = list(analyte = "clopidogrel", units = "mg", specimen = "tissue", verified = TRUE),
    wall_ileum = list(analyte = "clopidogrel", units = "mg", specimen = "tissue", verified = TRUE),
    wall_cecum = list(analyte = "clopidogrel", units = "mg", specimen = "tissue", verified = TRUE),
    wall_colon = list(analyte = "clopidogrel", units = "mg", specimen = "tissue", verified = TRUE),
    liver = list(analyte = "clopidogrel", units = "mg", specimen = "tissue", verified = TRUE),
    spleen = list(analyte = "clopidogrel", units = "mg", specimen = "tissue", verified = TRUE),
    lung = list(analyte = "clopidogrel", units = "mg", specimen = "tissue", verified = TRUE),
    heart = list(analyte = "clopidogrel", units = "mg", specimen = "tissue", verified = TRUE),
    brain = list(analyte = "clopidogrel", units = "mg", specimen = "tissue", verified = TRUE),
    muscle = list(analyte = "clopidogrel", units = "mg", specimen = "tissue", verified = TRUE),
    adipose = list(analyte = "clopidogrel", units = "mg", specimen = "tissue", verified = TRUE),
    skin = list(analyte = "clopidogrel", units = "mg", specimen = "tissue", verified = TRUE),
    kidney = list(analyte = "clopidogrel", units = "mg", specimen = "tissue", verified = TRUE),
    other = list(analyte = "clopidogrel", units = "mg", specimen = "tissue", verified = TRUE),
    venous = list(analyte = "clopidogrel", units = "mg", specimen = "whole blood", verified = TRUE),
    arterial = list(analyte = "clopidogrel", units = "mg", specimen = "whole blood", verified = TRUE),
    wall_stomach_oxoclop = list(analyte = "2-oxo-clopidogrel", units = "mg", specimen = "tissue", verified = TRUE),
    wall_duodenum_oxoclop = list(analyte = "2-oxo-clopidogrel", units = "mg", specimen = "tissue", verified = TRUE),
    wall_jejunum_oxoclop = list(analyte = "2-oxo-clopidogrel", units = "mg", specimen = "tissue", verified = TRUE),
    wall_ileum_oxoclop = list(analyte = "2-oxo-clopidogrel", units = "mg", specimen = "tissue", verified = TRUE),
    wall_cecum_oxoclop = list(analyte = "2-oxo-clopidogrel", units = "mg", specimen = "tissue", verified = TRUE),
    wall_colon_oxoclop = list(analyte = "2-oxo-clopidogrel", units = "mg", specimen = "tissue", verified = TRUE),
    liver_oxoclop = list(analyte = "2-oxo-clopidogrel", units = "mg", specimen = "tissue", verified = TRUE),
    spleen_oxoclop = list(analyte = "2-oxo-clopidogrel", units = "mg", specimen = "tissue", verified = TRUE),
    lung_oxoclop = list(analyte = "2-oxo-clopidogrel", units = "mg", specimen = "tissue", verified = TRUE),
    heart_oxoclop = list(analyte = "2-oxo-clopidogrel", units = "mg", specimen = "tissue", verified = TRUE),
    brain_oxoclop = list(analyte = "2-oxo-clopidogrel", units = "mg", specimen = "tissue", verified = TRUE),
    muscle_oxoclop = list(analyte = "2-oxo-clopidogrel", units = "mg", specimen = "tissue", verified = TRUE),
    adipose_oxoclop = list(analyte = "2-oxo-clopidogrel", units = "mg", specimen = "tissue", verified = TRUE),
    skin_oxoclop = list(analyte = "2-oxo-clopidogrel", units = "mg", specimen = "tissue", verified = TRUE),
    kidney_oxoclop = list(analyte = "2-oxo-clopidogrel", units = "mg", specimen = "tissue", verified = TRUE),
    other_oxoclop = list(analyte = "2-oxo-clopidogrel", units = "mg", specimen = "tissue", verified = TRUE),
    venous_oxoclop = list(analyte = "2-oxo-clopidogrel", units = "mg", specimen = "whole blood", verified = TRUE),
    arterial_oxoclop = list(analyte = "2-oxo-clopidogrel", units = "mg", specimen = "whole blood", verified = TRUE),
    wall_stomach_h4 = list(
      analyte = "clopidogrel active metabolite (H4 thiol)",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    ),
    wall_duodenum_h4 = list(
      analyte = "clopidogrel active metabolite (H4 thiol)",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    ),
    wall_jejunum_h4 = list(
      analyte = "clopidogrel active metabolite (H4 thiol)",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    ),
    wall_ileum_h4 = list(
      analyte = "clopidogrel active metabolite (H4 thiol)",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    ),
    wall_cecum_h4 = list(
      analyte = "clopidogrel active metabolite (H4 thiol)",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    ),
    wall_colon_h4 = list(
      analyte = "clopidogrel active metabolite (H4 thiol)",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    ),
    liver_h4 = list(
      analyte = "clopidogrel active metabolite (H4 thiol)",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    ),
    spleen_h4 = list(
      analyte = "clopidogrel active metabolite (H4 thiol)",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    ),
    lung_h4 = list(
      analyte = "clopidogrel active metabolite (H4 thiol)",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    ),
    heart_h4 = list(
      analyte = "clopidogrel active metabolite (H4 thiol)",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    ),
    brain_h4 = list(
      analyte = "clopidogrel active metabolite (H4 thiol)",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    ),
    muscle_h4 = list(
      analyte = "clopidogrel active metabolite (H4 thiol)",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    ),
    adipose_h4 = list(
      analyte = "clopidogrel active metabolite (H4 thiol)",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    ),
    skin_h4 = list(
      analyte = "clopidogrel active metabolite (H4 thiol)",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    ),
    kidney_h4 = list(
      analyte = "clopidogrel active metabolite (H4 thiol)",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    ),
    other_h4 = list(
      analyte = "clopidogrel active metabolite (H4 thiol)",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    ),
    venous_h4 = list(
      analyte = "clopidogrel active metabolite (H4 thiol)",
      units = "mg",
      specimen = "whole blood",
      verified = TRUE
    ),
    arterial_h4 = list(
      analyte = "clopidogrel active metabolite (H4 thiol)",
      units = "mg",
      specimen = "whole blood",
      verified = TRUE
    ),
    aggregation = list(
      analyte = "normalized maximal platelet aggregation",
      units = "fraction of baseline",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = NA_integer_,
    n_studies = NA_integer_,
    age_range = "adults (healthy volunteers and CAD patients in the cited clinical studies)",
    weight_range = "70 kg reference adult used for the physiology (Table 1 volumes sum to 70 L)",
    disease_state = "healthy volunteers; coronary artery disease with or without diabetes mellitus",
    dose_range = paste(
      "Oral clopidogrel 75-900 mg single doses and 300/75 mg or 600/150 mg",
      "loading/maintenance regimens (Tables 4 and 5)"
    ),
    regions = "Literature data from East Asia, Europe and North America",
    notes = paste(
      "Bottom-up prediction model, not a population fit to individual",
      "data. Structural parameters come from in vitro kinetics and the",
      "literature; kirre was fitted to one published IPA profile (Zhu",
      "2008). Predictions were compared against 15 healthy-volunteer PK",
      "reports, 12 healthy IPA reports, 5 CAD reports and 5 CAD + DM",
      "reports. The paper's visual predictive check estimated variances",
      "on Vmax,CYP2C19, CLint,CES1, Kt,i and kirre but did not report them,",
      "so the model is typical-value only."
    )
  )

  ini({
    # ---- Physiology: 70-kg adult (Table 1, 'Health' columns) -------------
    # Volumes in L, blood flows in L/min as printed; flows are converted to
    # L/h in model().
    v_lung <- fixed(0.94);     label("Lung volume (L)")               # Table 1 Lung
    v_heart <- fixed(0.27);    label("Heart volume (L)")              # Table 1 Heart
    v_brain <- fixed(1.53);    label("Brain volume (L)")              # Table 1 Brain
    v_muscle <- fixed(17.51);  label("Muscle volume (L)")             # Table 1 Muscle
    v_adipose <- fixed(22.20); label("Adipose volume (L)")            # Table 1 Adipose
    v_skin <- fixed(1.65);     label("Skin volume (L)")               # Table 1 Skin
    v_kidney <- fixed(0.23);   label("Kidney volume (L)")             # Table 1 Kidney
    v_spleen <- fixed(0.16);   label("Spleen volume (L)")             # Table 1 Spleen
    v_liver <- fixed(1.38);    label("Liver volume (L)")              # Table 1 Liver
    v_other <- fixed(17.75);   label("Rest-of-body volume (L)")       # Table 1 ROB
    v_venous <- fixed(1.91);   label("Venous blood volume (L)")       # Table 1 Vein
    v_arterial <- fixed(3.83); label("Arterial blood volume (L)")     # Table 1 Artery
    v_wall_stomach <- fixed(0.15);  label("Stomach wall volume (L)")  # Table 1 Stomach
    v_wall_duodenum <- fixed(0.02); label("Duodenum wall volume (L)") # Table 1 Duodenum
    v_wall_jejunum <- fixed(0.06);  label("Jejunum wall volume (L)")  # Table 1 Jejunum
    v_wall_ileum <- fixed(0.04);    label("Ileum wall volume (L)")    # Table 1 Ileum
    v_wall_cecum <- fixed(0.04);    label("Cecum wall volume (L)")    # Table 1 Caecum
    v_wall_colon <- fixed(0.34);    label("Colon wall volume (L)")    # Table 1 Colon

    q_co_min <- fixed(5.27);      label("Cardiac output = lung blood flow (L/min)")  # Table 1 Lung / Vein / Artery
    q_heart_min <- fixed(0.20);   label("Heart blood flow (L/min)")                  # Table 1 Heart
    q_brain_min <- fixed(0.80);   label("Brain blood flow (L/min)")                  # Table 1 Brain
    q_muscle_min <- fixed(0.54);  label("Muscle blood flow (L/min)")                 # Table 1 Muscle
    q_adipose_min <- fixed(0.46); label("Adipose blood flow (L/min)")                # Table 1 Adipose
    q_skin_min <- fixed(0.21);    label("Skin blood flow (L/min)")                   # Table 1 Skin
    q_kidney_min <- fixed(0.89);  label("Kidney blood flow (L/min)")                 # Table 1 Kidney
    q_spleen_min <- fixed(0.16);  label("Spleen blood flow (L/min)")                 # Table 1 Spleen
    q_ha_min <- fixed(0.42);      label("Hepatic-artery blood flow Qliv (L/min)")    # Table 1 Liver (Qliv of Eq 11)
    q_other_min <- fixed(0.68);   label("Rest-of-body blood flow (L/min)")           # Table 1 ROB
    q_wall_stomach_min <- fixed(0.13);  label("Stomach wall blood flow (L/min)")     # Table 1 Stomach
    q_wall_duodenum_min <- fixed(0.08); label("Duodenum wall blood flow (L/min)")    # Table 1 Duodenum
    q_wall_jejunum_min <- fixed(0.30);  label("Jejunum wall blood flow (L/min)")     # Table 1 Jejunum
    q_wall_ileum_min <- fixed(0.17);    label("Ileum wall blood flow (L/min)")       # Table 1 Ileum
    q_wall_cecum_min <- fixed(0.03);    label("Cecum wall blood flow (L/min)")       # Table 1 Caecum
    q_wall_colon_min <- fixed(0.20);    label("Colon wall blood flow (L/min)")       # Table 1 Colon

    # ---- Clopidogrel physicochemistry and Kt/p (Table 2, 'CLOP') ----------
    fu <- fixed(0.02);     label("Clopidogrel fraction unbound in plasma (unitless)")          # Table 2 fup
    bp <- fixed(0.57);     label("Clopidogrel blood-to-plasma ratio Rbp (unitless)")           # Table 2 Rbp
    fumic <- fixed(0.015); label("Clopidogrel fraction unbound in hepatic microsomes (unitless)") # Table 2 fumic
    fu_gut <- fixed(0.02); label("Clopidogrel fraction unbound in gut wall fugut (unitless)")   # Methods 'In gut lumen': fugut 0.02
    lkp_adipose <- fixed(log(6.634)); label("Clopidogrel adipose:plasma partition coefficient (unitless)") # Table 2 Adipose
    lkp_brain <- fixed(log(6.675));   label("Clopidogrel brain:plasma partition coefficient (unitless)")   # Table 2 Brain
    lkp_gut <- fixed(log(6.586));     label("Clopidogrel gut:plasma partition coefficient (unitless)")     # Table 2 Gut
    lkp_heart <- fixed(log(6.065));   label("Clopidogrel heart:plasma partition coefficient (unitless)")   # Table 2 Heart
    lkp_kidney <- fixed(log(5.309));  label("Clopidogrel kidney:plasma partition coefficient (unitless)")  # Table 2 Kidney
    lkp_liver <- fixed(log(5.625));   label("Clopidogrel liver:plasma partition coefficient (unitless)")   # Table 2 Liver
    lkp_lung <- fixed(log(3.566));    label("Clopidogrel lung:plasma partition coefficient (unitless)")    # Table 2 Lung
    lkp_muscle <- fixed(log(6.171));  label("Clopidogrel muscle:plasma partition coefficient (unitless)")  # Table 2 Muscle
    lkp_skin <- fixed(log(5.309));    label("Clopidogrel skin:plasma partition coefficient (unitless)")    # Table 2 Skin
    lkp_spleen <- fixed(log(5.182));  label("Clopidogrel spleen:plasma partition coefficient (unitless)")  # Table 2 Spleen
    lkp_other <- fixed(log(0.001));   label("Clopidogrel rest-of-body:plasma partition coefficient (unitless)") # Table 2 ROB (assumed)

    # ---- 2-oxo-clopidogrel (Table 2, '2-oxo-CLOP') ------------------------
    fu_oxoclop <- fixed(0.0742);   label("2-oxo-clopidogrel fraction unbound in plasma (unitless)")    # Table 2 fup
    bp_oxoclop <- fixed(0.68);     label("2-oxo-clopidogrel blood-to-plasma ratio (unitless)")         # Table 2 Rbp
    fumic_oxoclop <- fixed(0.180); label("2-oxo-clopidogrel fraction unbound in microsomes (unitless)") # Table 2 fumic
    lkp_adipose_oxoclop <- fixed(log(0.151)); label("2-oxo-clopidogrel adipose Kt/p (unitless)") # Table 2 Adipose
    lkp_brain_oxoclop <- fixed(log(8.784));   label("2-oxo-clopidogrel brain Kt/p (unitless)")   # Table 2 Brain
    lkp_gut_oxoclop <- fixed(log(4.866));     label("2-oxo-clopidogrel gut Kt/p (unitless)")     # Table 2 Gut
    lkp_heart_oxoclop <- fixed(log(6.864));   label("2-oxo-clopidogrel heart Kt/p (unitless)")   # Table 2 Heart
    lkp_kidney_oxoclop <- fixed(log(7.933));  label("2-oxo-clopidogrel kidney Kt/p (unitless)")  # Table 2 Kidney
    lkp_liver_oxoclop <- fixed(log(7.988));   label("2-oxo-clopidogrel liver Kt/p (unitless)")   # Table 2 Liver
    lkp_lung_oxoclop <- fixed(log(6.810));    label("2-oxo-clopidogrel lung Kt/p (unitless)")    # Table 2 Lung
    lkp_muscle_oxoclop <- fixed(log(9.525));  label("2-oxo-clopidogrel muscle Kt/p (unitless)")  # Table 2 Muscle
    lkp_skin_oxoclop <- fixed(log(1.390));    label("2-oxo-clopidogrel skin Kt/p (unitless)")    # Table 2 Skin
    lkp_spleen_oxoclop <- fixed(log(8.748));  label("2-oxo-clopidogrel spleen Kt/p (unitless)")  # Table 2 Spleen
    lkp_other_oxoclop <- fixed(log(0.001));   label("2-oxo-clopidogrel rest-of-body Kt/p (unitless)") # Table 2 ROB (assumed)

    # ---- CLOP-AM, the H4 active thiol (Table 2, 'CLOP-AM') ----------------
    fu_h4 <- fixed(0.0791); label("CLOP-AM fraction unbound in plasma (unitless)") # Table 2 fup
    bp_h4 <- fixed(0.58);   label("CLOP-AM blood-to-plasma ratio (unitless)")      # Table 2 Rbp
    lkp_adipose_h4 <- fixed(log(0.098)); label("CLOP-AM adipose Kt/p (unitless)") # Table 2 Adipose
    lkp_brain_h4 <- fixed(log(3.544));   label("CLOP-AM brain Kt/p (unitless)")   # Table 2 Brain
    lkp_gut_h4 <- fixed(log(1.811));     label("CLOP-AM gut Kt/p (unitless)")     # Table 2 Gut
    lkp_heart_h4 <- fixed(log(2.805));   label("CLOP-AM heart Kt/p (unitless)")   # Table 2 Heart
    lkp_kidney_h4 <- fixed(log(3.597));  label("CLOP-AM kidney Kt/p (unitless)")  # Table 2 Kidney
    lkp_liver_h4 <- fixed(log(3.640));   label("CLOP-AM liver Kt/p (unitless)")   # Table 2 Liver
    lkp_lung_h4 <- fixed(log(1.624));    label("CLOP-AM lung Kt/p (unitless)")    # Table 2 Lung
    lkp_muscle_h4 <- fixed(log(2.802));  label("CLOP-AM muscle Kt/p (unitless)")  # Table 2 Muscle
    lkp_skin_h4 <- fixed(log(0.601));    label("CLOP-AM skin Kt/p (unitless)")    # Table 2 Skin
    lkp_spleen_h4 <- fixed(log(3.243));  label("CLOP-AM spleen Kt/p (unitless)")  # Table 2 Spleen
    lkp_other_h4 <- fixed(log(0.001));   label("CLOP-AM rest-of-body Kt/p (unitless)") # Table 2 ROB (assumed)

    # ---- Molecular weights (not printed by Xu 2020) ------------------------
    # Needed to convert mass concentrations to the molar units of Km (uM)
    # and kirre (mL/nmol/h); doses are mg and outputs ng/mL.
    mw <- fixed(321.82);         label("Clopidogrel molecular weight (g/mol)")       # not in Xu 2020; Jung 2024 Methods 'PK model' (doi:10.1002/psp4.13053) 321.82
    mw_oxoclop <- fixed(337.82); label("2-oxo-clopidogrel molecular weight (g/mol)") # not in Xu 2020; C16H16ClNO3S = clopidogrel + one O (321.82 + 16.00)
    mw_h4 <- fixed(355.83);      label("CLOP-AM (H4) molecular weight (g/mol)")      # not in Xu 2020; Jung 2024 Methods 'PK model' 355.83

    # ---- Gastrointestinal transit and absorption (Methods, Eqs 2-8) -------
    ktr_stomach <- fixed(4.8);   label("Gastric emptying rate constant Kt,0, healthy/CAD (1/h)") # Methods 'In stomach': Kt,0 4.8 h-1
    ktr_duodenum <- fixed(4.2);  label("Duodenum transit rate constant Kt,1, healthy/CAD (1/h)") # Methods 'In gut lumen': 4.2
    ktr_jejunum <- fixed(1.8);   label("Jejunum transit rate constant Kt,2, healthy/CAD (1/h)")  # Methods 'In gut lumen': 1.8
    ktr_ileum <- fixed(2.4);     label("Ileum transit rate constant Kt,3, healthy/CAD (1/h)")    # Methods 'In gut lumen': 2.4
    ktr_cecum <- fixed(0.18);    label("Cecum transit rate constant Kt,4, healthy/CAD (1/h)")    # Methods 'In gut lumen': 0.18
    ktr_colon <- fixed(0.06);    label("Colon transit rate constant Kt,5, healthy/CAD (1/h)")    # Methods 'In gut lumen': 0.06
    ktr_stomach_dm <- fixed(2.31);  label("Gastric emptying rate constant in DM (1/h)")  # Methods 'CAD Patients With DM': 2.31
    ktr_duodenum_dm <- fixed(2.30); label("Duodenum transit rate constant in DM (1/h)")  # Methods 'CAD Patients With DM': 2.30
    ktr_jejunum_dm <- fixed(0.99);  label("Jejunum transit rate constant in DM (1/h)")   # Methods 'CAD Patients With DM': 0.99
    ktr_ileum_dm <- fixed(1.32);    label("Ileum transit rate constant in DM (1/h)")     # Methods 'CAD Patients With DM': 1.32
    ktr_cecum_dm <- fixed(0.20);    label("Cecum transit rate constant in DM (1/h)")     # Methods 'CAD Patients With DM': 0.20
    ktr_colon_dm <- fixed(0.04);    label("Colon transit rate constant in DM (1/h)")     # Methods 'CAD Patients With DM': 0.04
    ka_duodenum <- fixed(0.21); label("Duodenum absorption rate constant Ka,1 (1/h)")      # Methods: 'calculated Ka,i ... 0.21, 0.26, and 0.29 h-1'
    ka_jejunum <- fixed(0.26);  label("Jejunum absorption rate constant Ka,2 (1/h)")       # Methods: Ka,2 0.26 h-1
    ka_ileum <- fixed(0.29);    label("Ileum absorption rate constant Ka,3 (1/h)")         # Methods: Ka,3 0.29 h-1
    kef_duodenum <- fixed(0.07); label("Duodenum wall P-gp efflux rate constant Kb,1 (1/h)") # Methods: 'calculated Kb,i ... 0.07, 0.12, and 0.16 h-1'
    kef_jejunum <- fixed(0.12);  label("Jejunum wall P-gp efflux rate constant Kb,2 (1/h)")  # Methods: Kb,2 0.12 h-1
    kef_ileum <- fixed(0.16);    label("Ileum wall P-gp efflux rate constant Kb,3 (1/h)")    # Methods: Kb,3 0.16 h-1

    # ---- Hepatic metabolism (Eqs 11-15, Table 3) ---------------------------
    pbsf <- fixed(55120); label("Total hepatic microsomal protein PBSF (mg)") # Methods 'In liver': 55,120 mg
    # Clopidogrel -> 2-oxo-clopidogrel. Vmax in pmol/min/pmol P450, Km in uM
    # (total microsomal; the unbound Km is Km * fumic), content in pmol P450/mg.
    vmax_cyp1a2 <- fixed(2.27);  label("CYP1A2 Vmax, clopidogrel (pmol/min/pmol P450)")  # Table 3 CYP1A2 CLOP Vmax
    km_cyp1a2 <- fixed(1.58);    label("CYP1A2 Km, clopidogrel (uM)")                   # Table 3 CYP1A2 CLOP Km
    vmax_cyp2b6 <- fixed(7.66);  label("CYP2B6 Vmax, clopidogrel (pmol/min/pmol P450)")  # Table 3 CYP2B6 CLOP Vmax
    km_cyp2b6 <- fixed(2.08);    label("CYP2B6 Km, clopidogrel (uM)")                   # Table 3 CYP2B6 CLOP Km
    vmax_cyp2c19 <- fixed(7.52); label("CYP2C19 Vmax in EM, clopidogrel (pmol/min/pmol P450)") # Table 3 CYP2C19(EM) CLOP Vmax
    km_cyp2c19 <- fixed(1.12);   label("CYP2C19 Km, clopidogrel (uM)")                  # Table 3 CYP2C19 CLOP Km
    # 2-oxo-clopidogrel -> CLOP-AM.
    vmax_cyp2b6_oxoclop <- fixed(2.48);  label("CYP2B6 Vmax, 2-oxo-clopidogrel (pmol/min/pmol P450)")  # Table 3 CYP2B6 2-oxo-CLOP Vmax
    km_cyp2b6_oxoclop <- fixed(1.62);    label("CYP2B6 Km, 2-oxo-clopidogrel (uM)")                   # Table 3 CYP2B6 2-oxo-CLOP Km
    vmax_cyp2c9_oxoclop <- fixed(0.855); label("CYP2C9 Vmax, 2-oxo-clopidogrel (pmol/min/pmol P450)")  # Table 3 CYP2C9 2-oxo-CLOP Vmax
    km_cyp2c9_oxoclop <- fixed(18.1);    label("CYP2C9 Km, 2-oxo-clopidogrel (uM)")                   # Table 3 CYP2C9 2-oxo-CLOP Km
    vmax_cyp3a4_oxoclop <- fixed(3.63);  label("CYP3A4 Vmax, 2-oxo-clopidogrel (pmol/min/pmol P450)")  # Table 3 CYP3A4 2-oxo-CLOP Vmax
    km_cyp3a4_oxoclop <- fixed(27.8);    label("CYP3A4 Km, 2-oxo-clopidogrel (uM)")                   # Table 3 CYP3A4 2-oxo-CLOP Km
    vmax_cyp2c19_oxoclop <- fixed(9.06); label("CYP2C19 Vmax in EM, 2-oxo-clopidogrel (pmol/min/pmol P450)") # Table 3 CYP2C19(EM) 2-oxo-CLOP Vmax
    km_cyp2c19_oxoclop <- fixed(12.1);   label("CYP2C19 Km, 2-oxo-clopidogrel (uM)")                  # Table 3 CYP2C19 2-oxo-CLOP Km
    abund_cyp1a2 <- fixed(52);   label("CYP1A2 microsomal content (pmol P450/mg protein)")  # Table 3 enzyme content CYP1A2
    abund_cyp2b6 <- fixed(11);   label("CYP2B6 microsomal content (pmol P450/mg protein)")  # Table 3 enzyme content CYP2B6
    abund_cyp2c9 <- fixed(73);   label("CYP2C9 microsomal content (pmol P450/mg protein)")  # Table 3 enzyme content CYP2C9
    abund_cyp3a4 <- fixed(155);  label("CYP3A4 microsomal content (pmol P450/mg protein)")  # Table 3 enzyme content CYP3A4
    abund_cyp2c19 <- fixed(14);  label("CYP2C19 microsomal content (pmol P450/mg protein)") # Table 3 enzyme content CYP2C19
    clint_ces1 <- fixed(276650);  label("CES1 intrinsic clearance of clopidogrel (L/h)")        # Methods 'In liver': 276,650 (printed '1/h'; L/h per the 85% share, see vignette)
    clint_ces1_oxoclop <- fixed(2200); label("CES1 intrinsic clearance of 2-oxo-clopidogrel (L/h)") # Methods 'In liver': 2,200 l/h
    clint_ces1_h4 <- fixed(529);  label("CES1 intrinsic clearance of CLOP-AM (L/h)")           # Methods 'In liver': 529 l/h

    # ---- CYP2C19 phenotype, CAD and DM effects -----------------------------
    e_cyp2c19_um <- fixed(1.58); label("CYP2C19 activity in UM relative to EM (fold)") # Table 3 footnote b / Methods: 1.58-fold
    e_cyp2c19_im <- fixed(0.5);  label("CYP2C19 activity in IM relative to EM (fold)") # Table 3 footnote c / Methods: 50%
    e_cyp2c19_pm <- fixed(0);    label("CYP2C19 activity in PM relative to EM (fold)") # Table 3 footnote c / Methods: 0%
    e_cad_q <- fixed(0.90);      label("Blood-flow ratio Qtotal,CAD / Qtotal,health (fold)") # Eq 18; Methods: ratio 0.90 (Rerych 1978)
    e_cad_kirre <- fixed(0.7);   label("kirre in CAD relative to healthy (fold)")          # Methods: 'kirre value in CAD patients was corrected to 0.7 times'
    e_dm_cyp1a2 <- fixed(1.23);  label("CYP1A2 activity in DM relative to healthy (fold)")  # Methods 'CAD Patients With DM': 1.23
    e_dm_cyp2b6 <- fixed(0.55);  label("CYP2B6 activity in DM relative to healthy (fold)")  # Methods 'CAD Patients With DM': 0.55
    e_dm_cyp2c9 <- fixed(1.26);  label("CYP2C9 activity in DM relative to healthy (fold)")  # Methods 'CAD Patients With DM': 1.26
    e_dm_cyp2c19 <- fixed(0.54); label("CYP2C19 activity in DM relative to healthy (fold)") # Methods 'CAD Patients With DM': 0.54
    e_dm_cyp3a4 <- fixed(0.62);  label("CYP3A4 activity in DM relative to healthy (fold)")  # Methods 'CAD Patients With DM': 0.62
    e_dm_ces1 <- fixed(1.27);    label("CES1 activity in DM relative to healthy (fold)")    # Methods 'CAD Patients With DM': 1.27

    # ---- PD: platelet aggregation turnover (Eqs 16-17) ---------------------
    kout <- fixed(0.007804);  label("Platelet disaggregation rate constant kout (1/h)")       # Methods 'PD Kinetics': 0.007804 h-1 (platelet t1/2 3.7 days)
    m0 <- fixed(1);           label("Baseline normalized maximal platelet aggregation M0 (unitless)") # Methods 'PD Kinetics': M0 equal to 1
    kirre <- fixed(47.576);   label("CLOP-AM irreversible antiplatelet rate constant kirre, healthy (mL/nmol/h)") # Methods 'PD Kinetics': 47.576 ml/nmol/h

    # ---- Residual error ----------------------------------------------------
    # Deterministic prediction model; no residual error was reported, so the
    # SDs are fixed at zero rather than invented.
    propSd <- fixed(0);    label("Proportional residual error, clopidogrel (fraction; not reported)") # not reported
    propSd_h4 <- fixed(0); label("Proportional residual error, CLOP-AM (fraction; not reported)")     # not reported
    addSd_IPA <- fixed(0); label("Additive residual error, IPA (%; not reported)")                    # not reported
  })

  model({
    # ================= Disease and phenotype scaling ======================
    # CAD scales every blood flow by Qtotal,CAD / Qtotal,health (Eq 18).
    fq <- 1 + (e_cad_q - 1) * DIS_IHD
    q_co <- 60 * q_co_min * fq
    q_heart <- 60 * q_heart_min * fq
    q_brain <- 60 * q_brain_min * fq
    q_muscle <- 60 * q_muscle_min * fq
    q_adipose <- 60 * q_adipose_min * fq
    q_skin <- 60 * q_skin_min * fq
    q_kidney <- 60 * q_kidney_min * fq
    q_spleen <- 60 * q_spleen_min * fq
    q_ha <- 60 * q_ha_min * fq
    q_other <- 60 * q_other_min * fq
    q_wall_stomach <- 60 * q_wall_stomach_min * fq
    q_wall_duodenum <- 60 * q_wall_duodenum_min * fq
    q_wall_jejunum <- 60 * q_wall_jejunum_min * fq
    q_wall_ileum <- 60 * q_wall_ileum_min * fq
    q_wall_cecum <- 60 * q_wall_cecum_min * fq
    q_wall_colon <- 60 * q_wall_colon_min * fq
    # Liver outflow (Eq 11): Qliv + Qsp + Qst + sum(Qgw,i).
    q_liver_out <- q_ha + q_spleen + q_wall_stomach + q_wall_duodenum +
      q_wall_jejunum + q_wall_ileum + q_wall_cecum + q_wall_colon

    # DM replaces the transit rate constants (Methods 'CAD Patients With DM').
    kt0 <- ktr_stomach * (1 - DIS_DIAB) + ktr_stomach_dm * DIS_DIAB
    kt1 <- ktr_duodenum * (1 - DIS_DIAB) + ktr_duodenum_dm * DIS_DIAB
    kt2 <- ktr_jejunum * (1 - DIS_DIAB) + ktr_jejunum_dm * DIS_DIAB
    kt3 <- ktr_ileum * (1 - DIS_DIAB) + ktr_ileum_dm * DIS_DIAB
    kt4 <- ktr_cecum * (1 - DIS_DIAB) + ktr_cecum_dm * DIS_DIAB
    kt5 <- ktr_colon * (1 - DIS_DIAB) + ktr_colon_dm * DIS_DIAB

    # CYP2C19 phenotype (Table 3) and DM enzyme-activity multipliers.
    f_cyp2c19 <- (1 + (e_cyp2c19_um - 1) * CYP2C19_UM + (e_cyp2c19_im - 1) * CYP2C19_IM +
      (e_cyp2c19_pm - 1) * CYP2C19_PM) * (1 + (e_dm_cyp2c19 - 1) * DIS_DIAB)
    f_cyp1a2 <- 1 + (e_dm_cyp1a2 - 1) * DIS_DIAB
    f_cyp2b6 <- 1 + (e_dm_cyp2b6 - 1) * DIS_DIAB
    f_cyp2c9 <- 1 + (e_dm_cyp2c9 - 1) * DIS_DIAB
    f_cyp3a4 <- 1 + (e_dm_cyp3a4 - 1) * DIS_DIAB
    f_ces1 <- 1 + (e_dm_ces1 - 1) * DIS_DIAB

    # CAD lowers kirre to 0.7-fold; the DM value equals the CAD value.
    kirre_i <- kirre * (1 + (e_cad_kirre - 1) * DIS_IHD)
    kin <- kout * m0

    # ================= Concentrations =====================================
    # All states are amounts in mg and concentrations in mg/L. The
    # venous-equivalent concentration leaving tissue t is C_t / Kt/b.
    # ---- clopidogrel: tissue-to-blood partition coefficients Kt/b = Kt/p / Rbp ----
    kb_lung <- exp(lkp_lung) / bp
    kb_heart <- exp(lkp_heart) / bp
    kb_brain <- exp(lkp_brain) / bp
    kb_muscle <- exp(lkp_muscle) / bp
    kb_adipose <- exp(lkp_adipose) / bp
    kb_skin <- exp(lkp_skin) / bp
    kb_kidney <- exp(lkp_kidney) / bp
    kb_spleen <- exp(lkp_spleen) / bp
    kb_liver <- exp(lkp_liver) / bp
    kb_other <- exp(lkp_other) / bp
    kb_gut <- exp(lkp_gut) / bp
    kb_stomach <- kb_gut # stomach-wall Kt/p not printed; the gut value is assumed
    fub <- fu / bp

    c_arterial <- arterial / v_arterial
    c_venous <- venous / v_venous
    cv_lung <- lung / v_lung / kb_lung
    cv_heart <- heart / v_heart / kb_heart
    cv_brain <- brain / v_brain / kb_brain
    cv_muscle <- muscle / v_muscle / kb_muscle
    cv_adipose <- adipose / v_adipose / kb_adipose
    cv_skin <- skin / v_skin / kb_skin
    cv_kidney <- kidney / v_kidney / kb_kidney
    cv_spleen <- spleen / v_spleen / kb_spleen
    cv_liver <- liver / v_liver / kb_liver
    cv_other <- other / v_other / kb_other
    cv_wall_stomach <- wall_stomach / v_wall_stomach / kb_stomach
    c_wall_duodenum <- wall_duodenum / v_wall_duodenum
    cv_wall_duodenum <- c_wall_duodenum / kb_gut
    c_wall_jejunum <- wall_jejunum / v_wall_jejunum
    cv_wall_jejunum <- c_wall_jejunum / kb_gut
    c_wall_ileum <- wall_ileum / v_wall_ileum
    cv_wall_ileum <- c_wall_ileum / kb_gut
    c_wall_cecum <- wall_cecum / v_wall_cecum
    cv_wall_cecum <- c_wall_cecum / kb_gut
    c_wall_colon <- wall_colon / v_wall_colon
    cv_wall_colon <- c_wall_colon / kb_gut
    cu_liver <- cv_liver * fub

    # ---- 2-oxo-clopidogrel: tissue-to-blood partition coefficients Kt/b = Kt/p / Rbp ----
    kb_lung_oxoclop <- exp(lkp_lung_oxoclop) / bp_oxoclop
    kb_heart_oxoclop <- exp(lkp_heart_oxoclop) / bp_oxoclop
    kb_brain_oxoclop <- exp(lkp_brain_oxoclop) / bp_oxoclop
    kb_muscle_oxoclop <- exp(lkp_muscle_oxoclop) / bp_oxoclop
    kb_adipose_oxoclop <- exp(lkp_adipose_oxoclop) / bp_oxoclop
    kb_skin_oxoclop <- exp(lkp_skin_oxoclop) / bp_oxoclop
    kb_kidney_oxoclop <- exp(lkp_kidney_oxoclop) / bp_oxoclop
    kb_spleen_oxoclop <- exp(lkp_spleen_oxoclop) / bp_oxoclop
    kb_liver_oxoclop <- exp(lkp_liver_oxoclop) / bp_oxoclop
    kb_other_oxoclop <- exp(lkp_other_oxoclop) / bp_oxoclop
    kb_gut_oxoclop <- exp(lkp_gut_oxoclop) / bp_oxoclop
    kb_stomach_oxoclop <- kb_gut_oxoclop # stomach-wall Kt/p not printed; the gut value is assumed
    fub_oxoclop <- fu_oxoclop / bp_oxoclop

    c_arterial_oxoclop <- arterial_oxoclop / v_arterial
    c_venous_oxoclop <- venous_oxoclop / v_venous
    cv_lung_oxoclop <- lung_oxoclop / v_lung / kb_lung_oxoclop
    cv_heart_oxoclop <- heart_oxoclop / v_heart / kb_heart_oxoclop
    cv_brain_oxoclop <- brain_oxoclop / v_brain / kb_brain_oxoclop
    cv_muscle_oxoclop <- muscle_oxoclop / v_muscle / kb_muscle_oxoclop
    cv_adipose_oxoclop <- adipose_oxoclop / v_adipose / kb_adipose_oxoclop
    cv_skin_oxoclop <- skin_oxoclop / v_skin / kb_skin_oxoclop
    cv_kidney_oxoclop <- kidney_oxoclop / v_kidney / kb_kidney_oxoclop
    cv_spleen_oxoclop <- spleen_oxoclop / v_spleen / kb_spleen_oxoclop
    cv_liver_oxoclop <- liver_oxoclop / v_liver / kb_liver_oxoclop
    cv_other_oxoclop <- other_oxoclop / v_other / kb_other_oxoclop
    cv_wall_stomach_oxoclop <- wall_stomach_oxoclop / v_wall_stomach / kb_stomach_oxoclop
    c_wall_duodenum_oxoclop <- wall_duodenum_oxoclop / v_wall_duodenum
    cv_wall_duodenum_oxoclop <- c_wall_duodenum_oxoclop / kb_gut_oxoclop
    c_wall_jejunum_oxoclop <- wall_jejunum_oxoclop / v_wall_jejunum
    cv_wall_jejunum_oxoclop <- c_wall_jejunum_oxoclop / kb_gut_oxoclop
    c_wall_ileum_oxoclop <- wall_ileum_oxoclop / v_wall_ileum
    cv_wall_ileum_oxoclop <- c_wall_ileum_oxoclop / kb_gut_oxoclop
    c_wall_cecum_oxoclop <- wall_cecum_oxoclop / v_wall_cecum
    cv_wall_cecum_oxoclop <- c_wall_cecum_oxoclop / kb_gut_oxoclop
    c_wall_colon_oxoclop <- wall_colon_oxoclop / v_wall_colon
    cv_wall_colon_oxoclop <- c_wall_colon_oxoclop / kb_gut_oxoclop
    cu_liver_oxoclop <- cv_liver_oxoclop * fub_oxoclop

    # ---- CLOP-AM (H4): tissue-to-blood partition coefficients Kt/b = Kt/p / Rbp ----
    kb_lung_h4 <- exp(lkp_lung_h4) / bp_h4
    kb_heart_h4 <- exp(lkp_heart_h4) / bp_h4
    kb_brain_h4 <- exp(lkp_brain_h4) / bp_h4
    kb_muscle_h4 <- exp(lkp_muscle_h4) / bp_h4
    kb_adipose_h4 <- exp(lkp_adipose_h4) / bp_h4
    kb_skin_h4 <- exp(lkp_skin_h4) / bp_h4
    kb_kidney_h4 <- exp(lkp_kidney_h4) / bp_h4
    kb_spleen_h4 <- exp(lkp_spleen_h4) / bp_h4
    kb_liver_h4 <- exp(lkp_liver_h4) / bp_h4
    kb_other_h4 <- exp(lkp_other_h4) / bp_h4
    kb_gut_h4 <- exp(lkp_gut_h4) / bp_h4
    kb_stomach_h4 <- kb_gut_h4 # stomach-wall Kt/p not printed; the gut value is assumed
    fub_h4 <- fu_h4 / bp_h4

    c_arterial_h4 <- arterial_h4 / v_arterial
    c_venous_h4 <- venous_h4 / v_venous
    cv_lung_h4 <- lung_h4 / v_lung / kb_lung_h4
    cv_heart_h4 <- heart_h4 / v_heart / kb_heart_h4
    cv_brain_h4 <- brain_h4 / v_brain / kb_brain_h4
    cv_muscle_h4 <- muscle_h4 / v_muscle / kb_muscle_h4
    cv_adipose_h4 <- adipose_h4 / v_adipose / kb_adipose_h4
    cv_skin_h4 <- skin_h4 / v_skin / kb_skin_h4
    cv_kidney_h4 <- kidney_h4 / v_kidney / kb_kidney_h4
    cv_spleen_h4 <- spleen_h4 / v_spleen / kb_spleen_h4
    cv_liver_h4 <- liver_h4 / v_liver / kb_liver_h4
    cv_other_h4 <- other_h4 / v_other / kb_other_h4
    cv_wall_stomach_h4 <- wall_stomach_h4 / v_wall_stomach / kb_stomach_h4
    c_wall_duodenum_h4 <- wall_duodenum_h4 / v_wall_duodenum
    cv_wall_duodenum_h4 <- c_wall_duodenum_h4 / kb_gut_h4
    c_wall_jejunum_h4 <- wall_jejunum_h4 / v_wall_jejunum
    cv_wall_jejunum_h4 <- c_wall_jejunum_h4 / kb_gut_h4
    c_wall_ileum_h4 <- wall_ileum_h4 / v_wall_ileum
    cv_wall_ileum_h4 <- c_wall_ileum_h4 / kb_gut_h4
    c_wall_cecum_h4 <- wall_cecum_h4 / v_wall_cecum
    cv_wall_cecum_h4 <- c_wall_cecum_h4 / kb_gut_h4
    c_wall_colon_h4 <- wall_colon_h4 / v_wall_colon
    cv_wall_colon_h4 <- c_wall_colon_h4 / kb_gut_h4
    cu_liver_h4 <- cv_liver_h4 * fub_h4

    # ================= Gut absorption and efflux (Eqs 3, 9) ===============
    absorb_duodenum <- ka_duodenum * duodenum
    absorb_jejunum <- ka_jejunum * jejunum
    absorb_ileum <- ka_ileum * ileum
    efflux_duodenum <- kef_duodenum * c_wall_duodenum * v_wall_duodenum * fu_gut
    efflux_jejunum <- kef_jejunum * c_wall_jejunum * v_wall_jejunum * fu_gut
    efflux_ileum <- kef_ileum * c_wall_ileum * v_wall_ileum * fu_gut

    # ================= Hepatic metabolism (Eqs 11-15) =====================
    # CLint,CYP450 = sum over isoforms of Vmax * content / (Km * fumic + Cu),
    # in uL/min/mg protein; x PBSF (mg) x 60e-6 gives L/h. Km is in uM, so
    # the unbound liver concentration enters in uM (mg/L x 1000 / MW).
    cu_liver_um <- cu_liver * 1000 / mw
    cu_liver_oxoclop_um <- cu_liver_oxoclop * 1000 / mw_oxoclop
    clint_cyp_clop <- pbsf * 60e-6 * (
      vmax_cyp1a2 * abund_cyp1a2 * f_cyp1a2 / (km_cyp1a2 * fumic + cu_liver_um) +
        vmax_cyp2b6 * abund_cyp2b6 * f_cyp2b6 / (km_cyp2b6 * fumic + cu_liver_um) +
        vmax_cyp2c19 * abund_cyp2c19 * f_cyp2c19 / (km_cyp2c19 * fumic + cu_liver_um)
    )
    clint_cyp_oxoclop <- pbsf * 60e-6 * (
      vmax_cyp2b6_oxoclop * abund_cyp2b6 * f_cyp2b6 / (km_cyp2b6_oxoclop * fumic_oxoclop + cu_liver_oxoclop_um) +
        vmax_cyp2c9_oxoclop * abund_cyp2c9 * f_cyp2c9 / (km_cyp2c9_oxoclop * fumic_oxoclop + cu_liver_oxoclop_um) +
        vmax_cyp3a4_oxoclop * abund_cyp3a4 * f_cyp3a4 / (km_cyp3a4_oxoclop * fumic_oxoclop + cu_liver_oxoclop_um) +
        vmax_cyp2c19_oxoclop * abund_cyp2c19 * f_cyp2c19 / (km_cyp2c19_oxoclop * fumic_oxoclop + cu_liver_oxoclop_um)
    )
    # Metabolite formation adds the metabolized amount 1:1 BY MASS, exactly
    # as Eqs 13 and 15 are printed (no molecular-weight ratio); see the
    # vignette for why this as-run form is kept.
    met_cyp_clop <- clint_cyp_clop * cu_liver
    met_ces1_clop <- clint_ces1 * f_ces1 * cu_liver
    met_cyp_oxoclop <- clint_cyp_oxoclop * cu_liver_oxoclop
    met_ces1_oxoclop <- clint_ces1_oxoclop * f_ces1 * cu_liver_oxoclop
    met_ces1_h4 <- clint_ces1_h4 * f_ces1 * cu_liver_h4

    # ================= Stomach and gut lumen, clopidogrel (Eqs 2-3) =======
    d/dt(stomach) <- -kt0 * stomach
    d/dt(duodenum) <- kt0 * stomach - kt1 * duodenum - ka_duodenum * duodenum +
      efflux_duodenum
    d/dt(jejunum) <- kt1 * duodenum - kt2 * jejunum - ka_jejunum * jejunum +
      efflux_jejunum
    d/dt(ileum) <- kt2 * jejunum - kt3 * ileum - ka_ileum * ileum +
      efflux_ileum
    d/dt(cecum) <- kt3 * ileum - kt4 * cecum
    d/dt(colon) <- kt4 * cecum - kt5 * colon

    # ================= Clopidogrel (Eqs 1, 9, 11) ===========================
    d/dt(wall_stomach) <- q_wall_stomach * (c_arterial - cv_wall_stomach)
    d/dt(wall_duodenum) <- q_wall_duodenum * (c_arterial - cv_wall_duodenum) +
      absorb_duodenum - efflux_duodenum
    d/dt(wall_jejunum) <- q_wall_jejunum * (c_arterial - cv_wall_jejunum) +
      absorb_jejunum - efflux_jejunum
    d/dt(wall_ileum) <- q_wall_ileum * (c_arterial - cv_wall_ileum) +
      absorb_ileum - efflux_ileum
    d/dt(wall_cecum) <- q_wall_cecum * (c_arterial - cv_wall_cecum)
    d/dt(wall_colon) <- q_wall_colon * (c_arterial - cv_wall_colon)
    d/dt(heart) <- q_heart * (c_arterial - cv_heart)
    d/dt(brain) <- q_brain * (c_arterial - cv_brain)
    d/dt(muscle) <- q_muscle * (c_arterial - cv_muscle)
    d/dt(adipose) <- q_adipose * (c_arterial - cv_adipose)
    d/dt(skin) <- q_skin * (c_arterial - cv_skin)
    d/dt(kidney) <- q_kidney * (c_arterial - cv_kidney)
    d/dt(other) <- q_other * (c_arterial - cv_other)
    d/dt(spleen) <- q_spleen * (c_arterial - cv_spleen)
    d/dt(liver) <- q_ha * c_arterial + q_wall_stomach * cv_wall_stomach +
      q_spleen * cv_spleen + q_wall_duodenum * cv_wall_duodenum +
      q_wall_jejunum * cv_wall_jejunum + q_wall_ileum * cv_wall_ileum +
      q_wall_cecum * cv_wall_cecum + q_wall_colon * cv_wall_colon -
      q_liver_out * cv_liver - met_cyp_clop - met_ces1_clop
    d/dt(lung) <- q_co * (c_venous - cv_lung)
    d/dt(arterial) <- q_co * (cv_lung - c_arterial)
    d/dt(venous) <- q_heart * cv_heart + q_brain * cv_brain +
      q_muscle * cv_muscle + q_adipose * cv_adipose + q_skin * cv_skin +
      q_kidney * cv_kidney + q_other * cv_other + q_liver_out * cv_liver -
      q_co * c_venous

    # ================= 2-oxo-clopidogrel (Eqs 1, 10, 13) =====================
    d/dt(wall_stomach_oxoclop) <- q_wall_stomach * (c_arterial_oxoclop - cv_wall_stomach_oxoclop)
    d/dt(wall_duodenum_oxoclop) <- q_wall_duodenum * (c_arterial_oxoclop - cv_wall_duodenum_oxoclop)
    d/dt(wall_jejunum_oxoclop) <- q_wall_jejunum * (c_arterial_oxoclop - cv_wall_jejunum_oxoclop)
    d/dt(wall_ileum_oxoclop) <- q_wall_ileum * (c_arterial_oxoclop - cv_wall_ileum_oxoclop)
    d/dt(wall_cecum_oxoclop) <- q_wall_cecum * (c_arterial_oxoclop - cv_wall_cecum_oxoclop)
    d/dt(wall_colon_oxoclop) <- q_wall_colon * (c_arterial_oxoclop - cv_wall_colon_oxoclop)
    d/dt(heart_oxoclop) <- q_heart * (c_arterial_oxoclop - cv_heart_oxoclop)
    d/dt(brain_oxoclop) <- q_brain * (c_arterial_oxoclop - cv_brain_oxoclop)
    d/dt(muscle_oxoclop) <- q_muscle * (c_arterial_oxoclop - cv_muscle_oxoclop)
    d/dt(adipose_oxoclop) <- q_adipose * (c_arterial_oxoclop - cv_adipose_oxoclop)
    d/dt(skin_oxoclop) <- q_skin * (c_arterial_oxoclop - cv_skin_oxoclop)
    d/dt(kidney_oxoclop) <- q_kidney * (c_arterial_oxoclop - cv_kidney_oxoclop)
    d/dt(other_oxoclop) <- q_other * (c_arterial_oxoclop - cv_other_oxoclop)
    d/dt(spleen_oxoclop) <- q_spleen * (c_arterial_oxoclop - cv_spleen_oxoclop)
    d/dt(liver_oxoclop) <- q_ha * c_arterial_oxoclop + q_wall_stomach * cv_wall_stomach_oxoclop +
      q_spleen * cv_spleen_oxoclop + q_wall_duodenum * cv_wall_duodenum_oxoclop +
      q_wall_jejunum * cv_wall_jejunum_oxoclop + q_wall_ileum * cv_wall_ileum_oxoclop +
      q_wall_cecum * cv_wall_cecum_oxoclop + q_wall_colon * cv_wall_colon_oxoclop -
      q_liver_out * cv_liver_oxoclop + met_cyp_clop - met_cyp_oxoclop - met_ces1_oxoclop
    d/dt(lung_oxoclop) <- q_co * (c_venous_oxoclop - cv_lung_oxoclop)
    d/dt(arterial_oxoclop) <- q_co * (cv_lung_oxoclop - c_arterial_oxoclop)
    d/dt(venous_oxoclop) <- q_heart * cv_heart_oxoclop + q_brain * cv_brain_oxoclop +
      q_muscle * cv_muscle_oxoclop + q_adipose * cv_adipose_oxoclop + q_skin * cv_skin_oxoclop +
      q_kidney * cv_kidney_oxoclop + q_other * cv_other_oxoclop + q_liver_out * cv_liver_oxoclop -
      q_co * c_venous_oxoclop

    # ================= CLOP-AM / H4 (Eqs 1, 10, 15) ==========================
    d/dt(wall_stomach_h4) <- q_wall_stomach * (c_arterial_h4 - cv_wall_stomach_h4)
    d/dt(wall_duodenum_h4) <- q_wall_duodenum * (c_arterial_h4 - cv_wall_duodenum_h4)
    d/dt(wall_jejunum_h4) <- q_wall_jejunum * (c_arterial_h4 - cv_wall_jejunum_h4)
    d/dt(wall_ileum_h4) <- q_wall_ileum * (c_arterial_h4 - cv_wall_ileum_h4)
    d/dt(wall_cecum_h4) <- q_wall_cecum * (c_arterial_h4 - cv_wall_cecum_h4)
    d/dt(wall_colon_h4) <- q_wall_colon * (c_arterial_h4 - cv_wall_colon_h4)
    d/dt(heart_h4) <- q_heart * (c_arterial_h4 - cv_heart_h4)
    d/dt(brain_h4) <- q_brain * (c_arterial_h4 - cv_brain_h4)
    d/dt(muscle_h4) <- q_muscle * (c_arterial_h4 - cv_muscle_h4)
    d/dt(adipose_h4) <- q_adipose * (c_arterial_h4 - cv_adipose_h4)
    d/dt(skin_h4) <- q_skin * (c_arterial_h4 - cv_skin_h4)
    d/dt(kidney_h4) <- q_kidney * (c_arterial_h4 - cv_kidney_h4)
    d/dt(other_h4) <- q_other * (c_arterial_h4 - cv_other_h4)
    d/dt(spleen_h4) <- q_spleen * (c_arterial_h4 - cv_spleen_h4)
    d/dt(liver_h4) <- q_ha * c_arterial_h4 + q_wall_stomach * cv_wall_stomach_h4 +
      q_spleen * cv_spleen_h4 + q_wall_duodenum * cv_wall_duodenum_h4 +
      q_wall_jejunum * cv_wall_jejunum_h4 + q_wall_ileum * cv_wall_ileum_h4 +
      q_wall_cecum * cv_wall_cecum_h4 + q_wall_colon * cv_wall_colon_h4 -
      q_liver_out * cv_liver_h4 + met_cyp_oxoclop - met_ces1_h4
    d/dt(lung_h4) <- q_co * (c_venous_h4 - cv_lung_h4)
    d/dt(arterial_h4) <- q_co * (cv_lung_h4 - c_arterial_h4)
    d/dt(venous_h4) <- q_heart * cv_heart_h4 + q_brain * cv_brain_h4 +
      q_muscle * cv_muscle_h4 + q_adipose * cv_adipose_h4 + q_skin * cv_skin_h4 +
      q_kidney * cv_kidney_h4 + q_other * cv_other_h4 + q_liver_out * cv_liver_h4 -
      q_co * c_venous_h4

    # ================= PD: platelet aggregation (Eqs 16-17) ===============
    # kirre is per nmol/mL, so venous CLOP-AM is converted from mg/L.
    c_venous_h4_nmol <- c_venous_h4 * 1000 / mw_h4
    d/dt(aggregation) <- kin - kout * aggregation -
      kirre_i * c_venous_h4_nmol * fub_h4 * aggregation
    aggregation(0) <- m0

    # ================= Outputs ============================================
    # Venous plasma concentrations in ng/mL (blood mg/L / Rbp x 1000).
    Cc <- 1000 * c_venous / bp
    Cc_oxoclop <- 1000 * c_venous_oxoclop / bp_oxoclop
    Cc_h4 <- 1000 * c_venous_h4 / bp_h4
    IPA <- (1 - aggregation) * 100

    Cc ~ prop(propSd)
    Cc_h4 ~ prop(propSd_h4)
    IPA ~ add(addSd_IPA)
  })
}
