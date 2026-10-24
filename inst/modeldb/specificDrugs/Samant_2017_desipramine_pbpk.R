Samant_2017_desipramine_pbpk <- function() {
  description <- paste(
    "PBPK (whole-body, perfusion-limited, GastroPlus v8.6). Intravenous",
    "desipramine in a human reference adult: 13 perfusion-limited tissues",
    "(lung, adipose, muscle, liver, spleen, heart, brain, kidney, skin,",
    "reproductive organs, red marrow, yellow marrow, rest of body) between",
    "arterial and venous blood pools, with the gut and spleen draining into",
    "the liver. The paper's atypical distribution (Vss about 15 L/kg from",
    "phospholipid binding and lysosomal trapping) is carried by the",
    "deposited tissue:plasma partition coefficients, which the authors",
    "computed with a raised blood:plasma ratio before reverting it to 1.03.",
    "Renal clearance is fu x GFR. Hepatic clearance is anchored to the",
    "paper's own simulated total plasma clearance (Supplementary Table) and",
    "converted to a well-stirred intrinsic clearance inside model(); the",
    "paper's CYP2D6 in vitro scaling and ontogeny are not reproduced. The",
    "oral ACAT arm (enterocyte trapping) is not included. Default values are",
    "the 34.4-year, 83.4 kg male used for the 12.5 mg IV model-building arm;",
    "the deposited physiology and Kp sets for other ages and sexes are",
    "substituted through ini() (see the vignette). Deterministic",
    "typical-value model with no between-subject or residual variability.",
    sep = " "
  )
  reference <- paste(
    "Samant TS, Lukacova V, Schmidt S. Development and Qualification of",
    "Physiologically Based Pharmacokinetic Models for Drugs With Atypical",
    "Distribution Behavior: A Desipramine Case Study. CPT Pharmacometrics",
    "Syst Pharmacol. 2017;6(5):315-321. doi:10.1002/psp4.12180.",
    "Physiology, partition coefficients and observed data are taken from the",
    "deposited GastroPlus project in Supporting Information file",
    "PSP4-6-315-s003.zip (physiology file HumAmeMal34.4YO_83.4kg.pbk;",
    "compound database Desipramine-reordered.mdb, table PBPKParams, record",
    "'IV Adult 12.5 mg Male'); clearance, Vss and plasma unbound fraction by",
    "age are from the Supplementary Table (PSP4-6-315-s001.docx).",
    sep = " "
  )
  vignette <- "Samant_2017_desipramine_pbpk"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # GastroPlus doses an IV bolus into the venous-return blood pool
  # (physiology file Settings InfusionKey = 4 = 'Venous Return').
  dosing <- c("venous")

  # GastroPlus tissue names with no register entry: ReproOrg, RedMarrow
  # and YellowMarrow are separate perfusion-limited tissues in the
  # deposited physiology, each with its own volume, flow and Kp.
  paper_specific_compartments <- c("reproductive_organs", "red_marrow", "yellow_marrow")

  compartmentData <- list(
    arterial = list(analyte = "desipramine", units = "mg", specimen = "whole blood", verified = TRUE),
    venous = list(analyte = "desipramine", units = "mg", specimen = "whole blood", verified = TRUE),
    lung = list(analyte = "desipramine", units = "mg", specimen = "tissue", verified = TRUE),
    adipose = list(analyte = "desipramine", units = "mg", specimen = "tissue", verified = TRUE),
    muscle = list(analyte = "desipramine", units = "mg", specimen = "tissue", verified = TRUE),
    liver = list(analyte = "desipramine", units = "mg", specimen = "tissue", verified = TRUE),
    spleen = list(analyte = "desipramine", units = "mg", specimen = "tissue", verified = TRUE),
    heart = list(analyte = "desipramine", units = "mg", specimen = "tissue", verified = TRUE),
    brain = list(analyte = "desipramine", units = "mg", specimen = "tissue", verified = TRUE),
    kidney = list(analyte = "desipramine", units = "mg", specimen = "tissue", verified = TRUE),
    skin = list(analyte = "desipramine", units = "mg", specimen = "tissue", verified = TRUE),
    reproductive_organs = list(analyte = "desipramine", units = "mg", specimen = "tissue", verified = TRUE),
    red_marrow = list(analyte = "desipramine", units = "mg", specimen = "tissue", verified = TRUE),
    yellow_marrow = list(analyte = "desipramine", units = "mg", specimen = "tissue", verified = TRUE),
    other = list(analyte = "desipramine", units = "mg", specimen = "tissue", verified = TRUE)
  )

  # Physiology is set per individual through ini(), not through data
  # columns: the deposited project carries one complete physiology file per
  # simulated subject rather than covariate equations.
  covariateData <- list()

  population <- list(
    species = "human",
    n_subjects = NA_integer_,
    n_studies = 4L,
    age_range = "Adults 22-42 years (simulated reference individuals); paediatric extrapolation 0-15 years",
    weight_range = "70-85.1 kg adults; 3.1-69.4 kg paediatric reference individuals",
    disease_state = "Healthy volunteers (literature mean profiles)",
    dose_range = "12.5 and 50 mg IV bolus (this model); 50 and 100 mg oral (paper, not included)",
    regions = "Not reported",
    notes = paste(
      "No individual data were fitted. The model was built against the",
      "published mean profile after 12.5 mg IV in normal volunteers",
      "(Ciraulo 1988, paper ref 22) and qualified against 50 mg IV",
      "(Brosen and Gram 1988, ref 23); the oral arms used Kurtz 1997 (50 mg,",
      "ref 24) and Brosen 1986 (100 mg, ref 25). Each arm is simulated in one",
      "GastroPlus PEAR reference individual whose physiology file is",
      "deposited: 83.4 kg 34.4-year male for 12.5 mg IV and 70.1 kg 30-year",
      "female for 50 mg IV. Paediatric reference individuals of both sexes",
      "at birth, 0.5, 1, 2.5, 5, 10 and 15 years are deposited as well.",
      sep = " "
    )
  )

  ini({
    # ===== Drug-specific parameters =====
    bp <- fixed(1.03); label("Blood-to-plasma concentration ratio used in simulation (unitless)") # Table 1 'B/P ratio 1.03'; Methods: reverted to 1.03 after the Kp values were computed; mdb Drugs.RatioBloodPlasma = 1.03
    fup <- fixed(0.178); label("Fraction unbound in plasma, adjusted (fraction)") # Supplementary Table, 'Fup' column, adults = 17.8 (GastroPlus adjusted fup; Table 1 unadjusted fu = 0.19)
    cl_kg <- fixed(508); label("Total plasma clearance per kg body weight (mL/h/kg)") # Supplementary Table, 34.5 years, Males, 'CL (mL/h/kg)' = 508

    # Tissue:plasma partition coefficients for the 12.5 mg IV male arm,
    # computed by GastroPlus (Lukacova method) with the optimized male
    # blood:plasma ratio of 2.4 (Methods). mdb PBPKParams.Kp, record
    # 'IV Adult 12.5 mg Male'.
    kp_lung <- fixed(51.289); label("Tissue:plasma partition coefficient, lung (unitless)") # mdb PBPKParams 'IV Adult 12.5 mg Male' / Lung
    kp_adipose <- fixed(5.95113); label("Tissue:plasma partition coefficient, adipose (unitless)") # mdb PBPKParams 'IV Adult 12.5 mg Male' / Adipose
    kp_muscle <- fixed(20.4146); label("Tissue:plasma partition coefficient, muscle (unitless)") # mdb PBPKParams 'IV Adult 12.5 mg Male' / Muscle
    kp_liver <- fixed(60.0074); label("Tissue:plasma partition coefficient, liver (unitless)") # mdb PBPKParams 'IV Adult 12.5 mg Male' / Liver
    kp_spleen <- fixed(41.9188); label("Tissue:plasma partition coefficient, spleen (unitless)") # mdb PBPKParams 'IV Adult 12.5 mg Male' / Spleen
    kp_heart <- fixed(29.7683); label("Tissue:plasma partition coefficient, heart (unitless)") # mdb PBPKParams 'IV Adult 12.5 mg Male' / Heart
    kp_brain <- fixed(5.92563); label("Tissue:plasma partition coefficient, brain (unitless)") # mdb PBPKParams 'IV Adult 12.5 mg Male' / Brain
    kp_kidney <- fixed(66.019); label("Tissue:plasma partition coefficient, kidney (unitless)") # mdb PBPKParams 'IV Adult 12.5 mg Male' / Kidney
    kp_skin <- fixed(17.6216); label("Tissue:plasma partition coefficient, skin (unitless)") # mdb PBPKParams 'IV Adult 12.5 mg Male' / Skin
    kp_reproductive_organs <- fixed(66.019); label("Tissue:plasma partition coefficient, reproductive organs (unitless)") # mdb PBPKParams 'IV Adult 12.5 mg Male' / ReproOrg
    kp_red_marrow <- fixed(9.35359); label("Tissue:plasma partition coefficient, red marrow (unitless)") # mdb PBPKParams 'IV Adult 12.5 mg Male' / RedMarrow
    kp_yellow_marrow <- fixed(5.95113); label("Tissue:plasma partition coefficient, yellow marrow (unitless)") # mdb PBPKParams 'IV Adult 12.5 mg Male' / YellowMarrow
    kp_other <- fixed(41.9189); label("Tissue:plasma partition coefficient, rest of body (unitless)") # mdb PBPKParams 'IV Adult 12.5 mg Male' / RestOfBody

    # ===== System-specific parameters: 34.4-year, 83.4 kg male =====
    # Physiology file HumAmeMal34.4YO_83.4kg.pbk, kept in its own units
    # (volumes mL, perfusions mL/s); model() converts to L and L/h.
    wt <- fixed(83.4); label("Body weight (kg)") # pbk Properties Weight = 83.4; ipd 'Body Wt [kg] = 83.4'
    gfr <- fixed(1.94245796736); label("Glomerular filtration rate (mL/s)") # pbk GlomerularFiltrationRate = 1.94245796736

    vol_arterial <- fixed(2133.86); label("Arterial blood volume (mL)") # pbk Arterial Supply Volume = 2133.86269039111
    vol_venous <- fixed(4267.73); label("Venous blood volume (mL)") # pbk Venous Return Volume = 4267.72538078221
    vol_lung <- fixed(1089.46); label("Lung volume (mL)") # pbk Lung Volume = 1089.45945945946
    vol_adipose <- fixed(29799.5); label("Adipose volume (mL)") # pbk Adipose Volume = 29799.4910655919
    vol_muscle <- fixed(20286.5); label("Muscle volume (mL)") # pbk Muscle Volume = 20286.4864864865
    vol_liver <- fixed(1596.62); label("Liver volume (mL)") # pbk Liver Volume = 1596.62162162162
    vol_spleen <- fixed(162.687); label("Spleen volume (mL)") # pbk Spleen Volume = 162.687279369819
    vol_heart <- fixed(317.459); label("Heart volume (mL)") # pbk Heart Volume = 317.459182363375
    vol_brain <- fixed(1650.49); label("Brain volume (mL)") # pbk Brain Volume = 1650.49220193621
    vol_kidney <- fixed(271.923); label("Kidney volume (mL)") # pbk Kidney Volume = 271.922531422703
    vol_skin <- fixed(1915.95); label("Skin volume (mL)") # pbk Skin Volume = 1915.94594594594
    vol_reproductive_organs <- fixed(31.3584); label("Reproductive organ volume (mL)") # pbk ReproOrg Volume = 31.3584
    vol_red_marrow <- fixed(1142.18); label("Red marrow volume (mL)") # pbk RedMarrow Volume = 1142.1810379776
    vol_yellow_marrow <- fixed(3174.87); label("Yellow marrow volume (mL)") # pbk YellowMarrow Volume = 3174.86760664214
    vol_other <- fixed(11903.4); label("Rest-of-body volume (mL)") # pbk RestOfBody Volume = 11903.3962983905

    perf_hepatic_artery <- fixed(8.21551); label("Hepatic arterial blood flow (mL/s)") # pbk Hepatic Artery TissuePerfusion = 8.21551243642963
    perf_gut <- fixed(13.7877); label("Gut (portal) blood flow (mL/s)") # pbk ACAT Gut TissuePerfusion = 13.787726274472
    perf_adipose <- fixed(10.0467); label("Adipose blood flow (mL/s)") # pbk Adipose TissuePerfusion = 10.0467277668951
    perf_muscle <- fixed(10.2695); label("Muscle blood flow (mL/s)") # pbk Muscle TissuePerfusion = 10.2694789533951
    perf_spleen <- fixed(2.74525); label("Spleen blood flow (mL/s)") # pbk Spleen TissuePerfusion = 2.74525429822786
    perf_heart <- fixed(3.9106); label("Heart blood flow (mL/s)") # pbk Heart TissuePerfusion = 3.91059609545379
    perf_brain <- fixed(14.2038); label("Brain blood flow (mL/s)") # pbk Brain TissuePerfusion = 14.2037808812606
    perf_kidney <- fixed(16.8854); label("Kidney blood flow (mL/s)") # pbk Kidney TissuePerfusion = 16.8853844124118
    perf_skin <- fixed(3.87958); label("Skin blood flow (mL/s)") # pbk Skin TissuePerfusion = 3.87958059534539
    perf_reproductive_organs <- fixed(0.11112); label("Reproductive organ blood flow (mL/s)") # pbk ReproOrg TissuePerfusion = 0.111120246303024
    perf_red_marrow <- fixed(5.78198); label("Red marrow blood flow (mL/s)") # pbk RedMarrow TissuePerfusion = 5.78197894603607
    perf_yellow_marrow <- fixed(1.60719); label("Yellow marrow blood flow (mL/s)") # pbk YellowMarrow TissuePerfusion = 1.60719029995088
    perf_other <- fixed(6.02577); label("Rest-of-body blood flow (mL/s)") # pbk RestOfBody TissuePerfusion = 6.02576879381872
  })

  model({
    # ===== Unit conversion: mL -> L, mL/s -> L/h =====
    v_arterial <- vol_arterial / 1000
    v_venous <- vol_venous / 1000
    v_lung <- vol_lung / 1000
    v_adipose <- vol_adipose / 1000
    v_muscle <- vol_muscle / 1000
    v_liver <- vol_liver / 1000
    v_spleen <- vol_spleen / 1000
    v_heart <- vol_heart / 1000
    v_brain <- vol_brain / 1000
    v_kidney <- vol_kidney / 1000
    v_skin <- vol_skin / 1000
    v_reproductive_organs <- vol_reproductive_organs / 1000
    v_red_marrow <- vol_red_marrow / 1000
    v_yellow_marrow <- vol_yellow_marrow / 1000
    v_other <- vol_other / 1000

    q_hepatic_artery <- perf_hepatic_artery * 3.6
    q_gut <- perf_gut * 3.6
    q_adipose <- perf_adipose * 3.6
    q_muscle <- perf_muscle * 3.6
    q_spleen <- perf_spleen * 3.6
    q_heart <- perf_heart * 3.6
    q_brain <- perf_brain * 3.6
    q_kidney <- perf_kidney * 3.6
    q_skin <- perf_skin * 3.6
    q_reproductive_organs <- perf_reproductive_organs * 3.6
    q_red_marrow <- perf_red_marrow * 3.6
    q_yellow_marrow <- perf_yellow_marrow * 3.6
    q_other <- perf_other * 3.6

    # The liver receives the hepatic artery plus the portal outflow of gut
    # and spleen (pbk topology: all three ConnectOut = 'Hepatic Junction');
    # the deposited liver TissuePerfusion equals this sum.
    q_liver <- q_hepatic_artery + q_gut + q_spleen
    # Cardiac output: sum of the arterial-side flows. Equals the deposited
    # lung / arterial / venous TissuePerfusion in every physiology file.
    q_co <- q_hepatic_artery + q_gut + q_spleen + q_adipose + q_muscle +
      q_heart + q_brain + q_kidney + q_skin + q_reproductive_organs +
      q_red_marrow + q_yellow_marrow + q_other

    # ===== Clearance =====
    # Renal: fu x GFR on the plasma leaving the kidney (Methods: 'renal
    # clearance was estimated as fu*GFR'). In blood-flow terms the
    # intrinsic renal clearance is fup * GFR / BP.
    clint_renal <- fup * gfr * 3.6 / bp
    cl_renal_plasma <- bp * q_kidney * clint_renal / (q_kidney + clint_renal)
    # Hepatic: the total plasma clearance is taken from the paper's
    # Supplementary Table (cl_kg, mL/h/kg) and the renal share removed; the
    # remaining blood clearance is inverted through the well-stirred liver
    # model to the intrinsic clearance acting on liver venous blood. This
    # needs cl_hepatic_blood < q_liver.
    cl_total_plasma <- cl_kg / 1000 * wt
    cl_hepatic_blood <- (cl_total_plasma - cl_renal_plasma) / bp
    clint_hepatic <- q_liver * cl_hepatic_blood / (q_liver - cl_hepatic_blood)

    # ===== Concentrations (mg/L) =====
    c_arterial <- arterial / v_arterial
    c_venous <- venous / v_venous
    # Blood leaving a perfusion-limited tissue: c_tissue / Kp is the plasma
    # concentration, times BP for whole blood.
    cb_lung <- lung / v_lung / kp_lung * bp
    cb_adipose <- adipose / v_adipose / kp_adipose * bp
    cb_muscle <- muscle / v_muscle / kp_muscle * bp
    cb_liver <- liver / v_liver / kp_liver * bp
    cb_spleen <- spleen / v_spleen / kp_spleen * bp
    cb_heart <- heart / v_heart / kp_heart * bp
    cb_brain <- brain / v_brain / kp_brain * bp
    cb_kidney <- kidney / v_kidney / kp_kidney * bp
    cb_skin <- skin / v_skin / kp_skin * bp
    cb_reproductive_organs <- reproductive_organs / v_reproductive_organs / kp_reproductive_organs * bp
    cb_red_marrow <- red_marrow / v_red_marrow / kp_red_marrow * bp
    cb_yellow_marrow <- yellow_marrow / v_yellow_marrow / kp_yellow_marrow * bp
    cb_other <- other / v_other / kp_other * bp

    # ===== Mass balance (mg/h) =====
    d/dt(venous) <- q_adipose * cb_adipose + q_muscle * cb_muscle +
      q_liver * cb_liver + q_heart * cb_heart + q_brain * cb_brain +
      q_kidney * cb_kidney + q_skin * cb_skin +
      q_reproductive_organs * cb_reproductive_organs +
      q_red_marrow * cb_red_marrow + q_yellow_marrow * cb_yellow_marrow +
      q_other * cb_other - q_co * c_venous
    d/dt(lung) <- q_co * (c_venous - cb_lung)
    d/dt(arterial) <- q_co * (cb_lung - c_arterial)
    d/dt(adipose) <- q_adipose * (c_arterial - cb_adipose)
    d/dt(muscle) <- q_muscle * (c_arterial - cb_muscle)
    # Gut is a flow-through: the deposited physiology carries no gut tissue
    # volume (its parameters live in the ACAT module), so arterial blood
    # reaches the liver unchanged through the portal vein.
    d/dt(liver) <- q_hepatic_artery * c_arterial + q_gut * c_arterial +
      q_spleen * cb_spleen - q_liver * cb_liver - clint_hepatic * cb_liver
    d/dt(spleen) <- q_spleen * (c_arterial - cb_spleen)
    d/dt(heart) <- q_heart * (c_arterial - cb_heart)
    d/dt(brain) <- q_brain * (c_arterial - cb_brain)
    d/dt(kidney) <- q_kidney * (c_arterial - cb_kidney) - clint_renal * cb_kidney
    d/dt(skin) <- q_skin * (c_arterial - cb_skin)
    d/dt(reproductive_organs) <- q_reproductive_organs * (c_arterial - cb_reproductive_organs)
    d/dt(red_marrow) <- q_red_marrow * (c_arterial - cb_red_marrow)
    d/dt(yellow_marrow) <- q_yellow_marrow * (c_arterial - cb_yellow_marrow)
    d/dt(other) <- q_other * (c_arterial - cb_other)

    # ===== Observation =====
    # Venous plasma concentration (GastroPlus 'Venous Return' plasma),
    # mg/L x 1000 = ng/mL.
    Cc <- c_venous / bp * 1000
  })
}
