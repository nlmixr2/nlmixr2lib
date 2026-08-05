Gaohua_2023_permeabilityLimited_pbpk <- function() {
  description <- paste(
    "PBPK (whole-body, permeability-limited, 55 ODEs; bespoke",
    "MATLAB/SimBiology 2022a model). Theoretical what-if investigation of",
    "how passive permeability, metabolism, active transporters and dosing",
    "route determine the dynamics of the tissue/plasma partition",
    "coefficient Kp(t) and the volume of distribution Vd(t). Twelve",
    "tissues (adipose, bone, brain, heart, kidney, muscle, skin, liver,",
    "pancreas, spleen, gut, lung) each carry four subcompartments -",
    "residual blood cells (_bc), residual plasma (_plasma),",
    "extracellular water (_ew) and intracellular water (_iw) - and the",
    "three blood compartments (venous, arterial, portal vein) carry the",
    "blood-cell and plasma subcompartments only, giving 12*4 + 3*2 = 54",
    "disposition states plus a gut-lumen depot. Passive permeation runs",
    "between adjacent subcompartments (Eqs 8-10), active uptake and efflux",
    "transporters sit on the cell membrane between _ew and _iw (Eqs",
    "11-12), and metabolism may occur in any subcompartment (Eqs 13-16).",
    "There is no drug: the compound is a generic small molecule with fu =",
    "fi = 1 everywhere, so a perfusion-limited model would give Kp = 1 in",
    "every tissue and Vdss = 1 L/kg exactly; any departure from those",
    "values measures the permeability-limited model's divergence from the",
    "perfusion-limited one. Every permeability, metabolic clearance and",
    "transporter clearance is parameterised as a fold-multiple of that",
    "tissue's plasma flow (the fold_* parameters), which is exactly how",
    "the paper drives its what-if scenarios - each scenario in Tables 1-4",
    "is a one-parameter change. Deterministic typical-value simulation",
    "model: the paper reports no IIV and no residual-error model."
  )
  reference <- paste(
    "Gaohua L, Zhang M, Sychterz C, Chang M, Schmidt BJ. The Interplay of",
    "Permeability, Metabolism, Transporters, and Dosing in Determining the",
    "Dynamics of the Tissue/Plasma Partition Coefficient and Volume of",
    "Distribution - A Theoretical Investigation Using",
    "Permeability-Limited, Physiologically Based Pharmacokinetic",
    "Modeling. Int J Mol Sci. 2023;24(22):16224.",
    "doi:10.3390/ijms242216224.",
    sep = " "
  )
  vignette <- "Gaohua_2023_permeabilityLimited_pbpk"
  units    <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # No covariates: the paper simulates a single 70 kg reference adult and
  # varies only drug-specific fold multipliers.
  covariateData <- list()

  compartmentData <- list(
    adipose_bc         = list(analyte = "generic small molecule", units = "mg", specimen = "blood cell", verified = TRUE),
    adipose_plasma     = list(analyte = "generic small molecule", units = "mg", specimen = "plasma", verified = TRUE),
    adipose_ew         = list(analyte = "generic small molecule", units = "mg", specimen = "tissue", verified = TRUE),
    adipose_iw         = list(analyte = "generic small molecule", units = "mg", specimen = "tissue", verified = TRUE),
    bone_bc            = list(analyte = "generic small molecule", units = "mg", specimen = "blood cell", verified = TRUE),
    bone_plasma        = list(analyte = "generic small molecule", units = "mg", specimen = "plasma", verified = TRUE),
    bone_ew            = list(analyte = "generic small molecule", units = "mg", specimen = "tissue", verified = TRUE),
    bone_iw            = list(analyte = "generic small molecule", units = "mg", specimen = "tissue", verified = TRUE),
    brain_bc           = list(analyte = "generic small molecule", units = "mg", specimen = "blood cell", verified = TRUE),
    brain_plasma       = list(analyte = "generic small molecule", units = "mg", specimen = "plasma", verified = TRUE),
    brain_ew           = list(analyte = "generic small molecule", units = "mg", specimen = "tissue", verified = TRUE),
    brain_iw           = list(analyte = "generic small molecule", units = "mg", specimen = "tissue", verified = TRUE),
    heart_bc           = list(analyte = "generic small molecule", units = "mg", specimen = "blood cell", verified = TRUE),
    heart_plasma       = list(analyte = "generic small molecule", units = "mg", specimen = "plasma", verified = TRUE),
    heart_ew           = list(analyte = "generic small molecule", units = "mg", specimen = "tissue", verified = TRUE),
    heart_iw           = list(analyte = "generic small molecule", units = "mg", specimen = "tissue", verified = TRUE),
    kidney_bc          = list(analyte = "generic small molecule", units = "mg", specimen = "blood cell", verified = TRUE),
    kidney_plasma      = list(analyte = "generic small molecule", units = "mg", specimen = "plasma", verified = TRUE),
    kidney_ew          = list(analyte = "generic small molecule", units = "mg", specimen = "tissue", verified = TRUE),
    kidney_iw          = list(analyte = "generic small molecule", units = "mg", specimen = "tissue", verified = TRUE),
    muscle_bc          = list(analyte = "generic small molecule", units = "mg", specimen = "blood cell", verified = TRUE),
    muscle_plasma      = list(analyte = "generic small molecule", units = "mg", specimen = "plasma", verified = TRUE),
    muscle_ew          = list(analyte = "generic small molecule", units = "mg", specimen = "tissue", verified = TRUE),
    muscle_iw          = list(analyte = "generic small molecule", units = "mg", specimen = "tissue", verified = TRUE),
    skin_bc            = list(analyte = "generic small molecule", units = "mg", specimen = "blood cell", verified = TRUE),
    skin_plasma        = list(analyte = "generic small molecule", units = "mg", specimen = "plasma", verified = TRUE),
    skin_ew            = list(analyte = "generic small molecule", units = "mg", specimen = "tissue", verified = TRUE),
    skin_iw            = list(analyte = "generic small molecule", units = "mg", specimen = "tissue", verified = TRUE),
    liver_bc           = list(analyte = "generic small molecule", units = "mg", specimen = "blood cell", verified = TRUE),
    liver_plasma       = list(analyte = "generic small molecule", units = "mg", specimen = "plasma", verified = TRUE),
    liver_ew           = list(analyte = "generic small molecule", units = "mg", specimen = "tissue", verified = TRUE),
    liver_iw           = list(analyte = "generic small molecule", units = "mg", specimen = "tissue", verified = TRUE),
    pancreas_bc        = list(analyte = "generic small molecule", units = "mg", specimen = "blood cell", verified = TRUE),
    pancreas_plasma    = list(analyte = "generic small molecule", units = "mg", specimen = "plasma", verified = TRUE),
    pancreas_ew        = list(analyte = "generic small molecule", units = "mg", specimen = "tissue", verified = TRUE),
    pancreas_iw        = list(analyte = "generic small molecule", units = "mg", specimen = "tissue", verified = TRUE),
    spleen_bc          = list(analyte = "generic small molecule", units = "mg", specimen = "blood cell", verified = TRUE),
    spleen_plasma      = list(analyte = "generic small molecule", units = "mg", specimen = "plasma", verified = TRUE),
    spleen_ew          = list(analyte = "generic small molecule", units = "mg", specimen = "tissue", verified = TRUE),
    spleen_iw          = list(analyte = "generic small molecule", units = "mg", specimen = "tissue", verified = TRUE),
    gut_bc             = list(analyte = "generic small molecule", units = "mg", specimen = "blood cell", verified = TRUE),
    gut_plasma         = list(analyte = "generic small molecule", units = "mg", specimen = "plasma", verified = TRUE),
    gut_ew             = list(analyte = "generic small molecule", units = "mg", specimen = "tissue", verified = TRUE),
    gut_iw             = list(analyte = "generic small molecule", units = "mg", specimen = "tissue", verified = TRUE),
    lung_bc            = list(analyte = "generic small molecule", units = "mg", specimen = "blood cell", verified = TRUE),
    lung_plasma        = list(analyte = "generic small molecule", units = "mg", specimen = "plasma", verified = TRUE),
    lung_ew            = list(analyte = "generic small molecule", units = "mg", specimen = "tissue", verified = TRUE),
    lung_iw            = list(analyte = "generic small molecule", units = "mg", specimen = "tissue", verified = TRUE),
    venous_bc          = list(analyte = "generic small molecule", units = "mg", specimen = "blood cell", verified = TRUE),
    venous_plasma      = list(analyte = "generic small molecule", units = "mg", specimen = "plasma", verified = TRUE),
    arterial_bc        = list(analyte = "generic small molecule", units = "mg", specimen = "blood cell", verified = TRUE),
    arterial_plasma    = list(analyte = "generic small molecule", units = "mg", specimen = "plasma", verified = TRUE),
    portal_bc          = list(analyte = "generic small molecule", units = "mg", specimen = "blood cell", verified = TRUE),
    portal_plasma      = list(analyte = "generic small molecule", units = "mg", specimen = "plasma", verified = TRUE),
    depot              = list(analyte = "generic small molecule", units = "mg", specimen = "administration site", verified = TRUE)
  )

  population <- list(
    species     = "human",
    n_subjects  = 0,
    notes       = paste(
      "No subjects: a theoretical modelling exercise, not a fit to data.",
      "System physiology is a single 70 kg reference adult with a 300 L/h",
      "cardiac output and haematocrit 0.45 (paper Section 4.3.1 and Table",
      "A1); tissue volumes and blood flows are drawn from the Simcyp",
      "Simulator V21, ICRP and Brown et al., with minor adjustments so",
      "that both columns of Table A1 balance to exactly 100%. The",
      "compound is a hypothetical generic small molecule (Section",
      "4.3.2).",
      sep = " "
    )
  )

  ini({
    # ================== System-specific parameters ==================
    # Paper Section 4.3.1 and Table A1. Density 1 kg/L is assumed for all
    # tissues and blood, so tissue volumes in L are numerically equal to
    # their share of body weight in kg.
    bw      <- fixed(70);   label("Body weight (kg)")                          # Section 4.3.1 / Table A1 header (70 kg)
    qc      <- fixed(300);  label("Cardiac output (L/h)")                      # Section 4.3.1 / Table A1 header (300 L/h)
    hct     <- fixed(0.45); label("Haematocrit (fraction of blood that is blood cells)")  # Section 4.3.1 (45% of residual blood is blood cells)
    density <- fixed(1);    label("Tissue and blood density (kg/L)")           # Section 4.3.1 (density of 1 kg/L assumed for all tissues and blood)

    # Tissue volume as a fraction of body weight (Table A1, column 1).
    fvol_adipose   <- fixed(0.25    ); label("Fractional volume, adipose (L/kg body weight)")  # Table A1 (25% of 70 kg *, adjusted to balance the 70 kg total)
    fvol_bone      <- fixed(0.15    ); label("Fractional volume, bone (L/kg body weight)")  # Table A1 (15% of 70 kg *, adjusted to balance the 70 kg total)
    fvol_brain     <- fixed(0.02    ); label("Fractional volume, brain (L/kg body weight)")  # Table A1 (2% of 70 kg)
    fvol_heart     <- fixed(0.005   ); label("Fractional volume, heart (L/kg body weight)")  # Table A1 (0.5% of 70 kg)
    fvol_kidney    <- fixed(0.0045  ); label("Fractional volume, kidney (L/kg body weight)")  # Table A1 (0.45% of 70 kg)
    fvol_muscle    <- fixed(0.4     ); label("Fractional volume, muscle (L/kg body weight)")  # Table A1 (40% of 70 kg)
    fvol_skin      <- fixed(0.0371  ); label("Fractional volume, skin (L/kg body weight)")  # Table A1 (3.71% of 70 kg)
    fvol_liver     <- fixed(0.0257  ); label("Fractional volume, liver (L/kg body weight)")  # Table A1 (2.57% of 70 kg)
    fvol_pancreas  <- fixed(0.0014  ); label("Fractional volume, pancreas (L/kg body weight)")  # Table A1 (0.14% of 70 kg)
    fvol_spleen    <- fixed(0.0026  ); label("Fractional volume, spleen (L/kg body weight)")  # Table A1 (0.26% of 70 kg)
    fvol_gut       <- fixed(0.0171  ); label("Fractional volume, gut (L/kg body weight)")  # Table A1 (1.71% of 70 kg)
    fvol_lung      <- fixed(0.0076  ); label("Fractional volume, lung (L/kg body weight)")  # Table A1 (0.76% of 70 kg)
    fvol_blood     <- fixed(0.079); label("Fractional volume, total blood (L/kg body weight)")  # Table A1 (blood 7.9% of 70 kg)
    fven_blood     <- fixed(0.6666667); label("Fraction of total blood volume held in the venous compartment (rest arterial)")  # Table A1 (1/3 in arterial and 2/3 in venous)
    v_portal       <- fixed(0.008);     label("Portal vein blood volume (L)")  # supplement BasePBPK.sbproj compartment PV = 0.008 L (outside the Table A1 70 kg balance)

    # Blood flow as a fraction of cardiac output (Table A1, column 2). The
    # portal-vein flow is derived (pancreas + spleen + gut = 19%) and the
    # lung receives the whole cardiac output, so neither is a parameter.
    fq_adipose     <- fixed(0.05   ); label("Fractional blood flow, adipose (fraction of cardiac output)")  # Table A1 (5%)
    fq_bone        <- fixed(0.05   ); label("Fractional blood flow, bone (fraction of cardiac output)")  # Table A1 (5%)
    fq_brain       <- fixed(0.12   ); label("Fractional blood flow, brain (fraction of cardiac output)")  # Table A1 (12%)
    fq_heart       <- fixed(0.04   ); label("Fractional blood flow, heart (fraction of cardiac output)")  # Table A1 (4%)
    fq_kidney      <- fixed(0.19   ); label("Fractional blood flow, kidney (fraction of cardiac output)")  # Table A1 (19%)
    fq_muscle      <- fixed(0.245  ); label("Fractional blood flow, muscle (fraction of cardiac output)")  # Table A1 (24.5% *, adjusted to balance the 300 L/h cardiac output)
    fq_skin        <- fixed(0.05   ); label("Fractional blood flow, skin (fraction of cardiac output)")  # Table A1 (5%)
    fq_pancreas    <- fixed(0.01   ); label("Fractional blood flow, pancreas (fraction of cardiac output)")  # Table A1 (1%)
    fq_spleen      <- fixed(0.03   ); label("Fractional blood flow, spleen (fraction of cardiac output)")  # Table A1 (3%)
    fq_gut         <- fixed(0.15   ); label("Fractional blood flow, gut (fraction of cardiac output)")  # Table A1 (15%)
    fq_liver_arterial <- fixed(0.065); label("Fractional blood flow, hepatic artery (fraction of cardiac output)")  # Table A1 (liver 25.5%: arterial 6.5%, portal vein 19%)

    # Residual blood volume as a fraction of tissue volume, from blood
    # remaining in rat tissues after bleeding (paper Section 4.3.1, ref
    # [33]); split into _bc and _plasma by haematocrit.
    frb_adipose    <- fixed(0.00625 ); label("Residual blood volume fraction, adipose (L/L tissue)")  # supplement BasePBPK.sbproj rule rbc_adipose = 0.45 * 0.00625 * Vadipose
    frb_bone       <- fixed(0.0191  ); label("Residual blood volume fraction, bone (L/L tissue)")  # supplement BasePBPK.sbproj rule rbc_bone = 0.45 * 0.0191 * Vbone
    frb_brain      <- fixed(0.0135  ); label("Residual blood volume fraction, brain (L/L tissue)")  # supplement BasePBPK.sbproj rule rbc_brain = 0.45 * 0.0135 * Vbrain
    frb_heart      <- fixed(0.061   ); label("Residual blood volume fraction, heart (L/L tissue)")  # supplement BasePBPK.sbproj rule rbc_heart = 0.45 * 0.061 * Vheart
    frb_kidney     <- fixed(0.0458  ); label("Residual blood volume fraction, kidney (L/L tissue)")  # supplement BasePBPK.sbproj rule rbc_kidney = 0.45 * 0.0458 * Vkidney
    frb_muscle     <- fixed(0.004   ); label("Residual blood volume fraction, muscle (L/L tissue)")  # supplement BasePBPK.sbproj rule rbc_muscle = 0.45 * 0.004 * Vmuscle
    frb_skin       <- fixed(0.0021  ); label("Residual blood volume fraction, skin (L/L tissue)")  # supplement BasePBPK.sbproj rule rbc_skin = 0.45 * 0.0021 * Vskin
    frb_liver      <- fixed(0.0572  ); label("Residual blood volume fraction, liver (L/L tissue)")  # supplement BasePBPK.sbproj rule rbc_liver = 0.45 * 0.0572 * Vliver
    frb_pancreas   <- fixed(0.0321  ); label("Residual blood volume fraction, pancreas (L/L tissue)")  # supplement BasePBPK.sbproj rule rbc_pancreas = 0.45 * 0.0321 * Vpancreas
    frb_spleen     <- fixed(0.321   ); label("Residual blood volume fraction, spleen (L/L tissue)")  # supplement BasePBPK.sbproj rule rbc_spleen = 0.45 * 0.321 * Vspleen
    frb_gut        <- fixed(0.0104  ); label("Residual blood volume fraction, gut (L/L tissue)")  # supplement BasePBPK.sbproj rule rbc_gut = 0.45 * 0.0104 * Vgut
    frb_lung       <- fixed(0.175   ); label("Residual blood volume fraction, lung (L/L tissue)")  # supplement BasePBPK.sbproj rule rbc_lung = 0.45 * 0.175 * Vlung

    # Extracellular water as a fraction of tissue volume, as implemented in
    # the Simcyp Simulator V21 (paper Section 4.1, ref [14]). Intracellular
    # water is the remainder of the tissue volume.
    few_adipose    <- fixed(0.141 ); label("Extracellular water volume fraction, adipose (L/L tissue)")  # supplement BasePBPK.sbproj rule ew_adipose = 0.141 * Vadipose
    few_bone       <- fixed(0.098 ); label("Extracellular water volume fraction, bone (L/L tissue)")  # supplement BasePBPK.sbproj rule ew_bone = 0.098 * Vbone
    few_brain      <- fixed(0.092 ); label("Extracellular water volume fraction, brain (L/L tissue)")  # supplement BasePBPK.sbproj rule ew_brain = 0.092 * Vbrain
    few_heart      <- fixed(0.313 ); label("Extracellular water volume fraction, heart (L/L tissue)")  # supplement BasePBPK.sbproj rule ew_heart = 0.313 * Vheart
    few_kidney     <- fixed(0.283 ); label("Extracellular water volume fraction, kidney (L/L tissue)")  # supplement BasePBPK.sbproj rule ew_kidney = 0.283 * Vkidney
    few_muscle     <- fixed(0.091 ); label("Extracellular water volume fraction, muscle (L/L tissue)")  # supplement BasePBPK.sbproj rule ew_muscle = 0.091 * Vmuscle
    few_skin       <- fixed(0.623 ); label("Extracellular water volume fraction, skin (L/L tissue)")  # supplement BasePBPK.sbproj rule ew_skin = 0.623 * Vskin
    few_liver      <- fixed(0.165 ); label("Extracellular water volume fraction, liver (L/L tissue)")  # supplement BasePBPK.sbproj rule ew_liver = 0.165 * Vliver
    few_pancreas   <- fixed(0.12  ); label("Extracellular water volume fraction, pancreas (L/L tissue)")  # supplement BasePBPK.sbproj rule ew_pancreas = 0.12 * Vpancreas
    few_spleen     <- fixed(0.208 ); label("Extracellular water volume fraction, spleen (L/L tissue)")  # supplement BasePBPK.sbproj rule ew_spleen = 0.208 * Vspleen
    few_gut        <- fixed(0.267 ); label("Extracellular water volume fraction, gut (L/L tissue)")  # supplement BasePBPK.sbproj rule ew_gut = 0.267 * Vgut
    few_lung       <- fixed(0.348 ); label("Extracellular water volume fraction, lung (L/L tissue)")  # supplement BasePBPK.sbproj rule ew_lung = 0.348 * Vlung

    # ================== Compound-specific parameters ==================
    # Paper Section 4.4.1: the residual plasma and blood cells are assumed
    # to be in immediate equilibrium, so their mutual PS is set very high.
    ps_bc_plasma <- fixed(1000); label("Passive permeability-surface area between residual blood cells and residual plasma (L/h)")  # Section 4.4.1 (PStc/tp were set to 1000 L/h)

    # The eight what-if multipliers. Every PS, metabolic clearance and
    # transporter clearance is expressed as a fold-multiple of that
    # tissue's plasma flow (blood-cell flow for the blood-cell metabolic
    # clearance), exactly as the paper drives its scenarios: "varying the
    # PS ... around the tissue blood flow" (Section 4.4.1), "varying the
    # metabolic clearance (CLmet) around tissue blood flow" (Section
    # 4.4.2), and "varying the uptake or efflux transporter clearance
    # ... around the tissue blood flow" (Section 4.4.3). The shipped
    # defaults are the base model as supplied: 1% of the tissue blood flow
    # for all eight (Section 4.5 and the supplement BasePBPK.sbproj).
    fold_ps_plasma_ew         <- fixed(0.01); label("PS between residual plasma and extracellular water, as a fold of tissue plasma flow")  # Section 4.5 / supplement BasePBPK.sbproj parameter FoldPlasmaEw = 0.01
    fold_ps_ew_iw             <- fixed(0.01); label("PS between extracellular and intracellular water, as a fold of tissue plasma flow")  # Section 4.5 / supplement BasePBPK.sbproj parameter FoldEwIw = 0.01
    fold_uptake               <- fixed(0.01); label("Active uptake transporter clearance (ew to iw), as a fold of tissue plasma flow")  # Section 4.5 / supplement BasePBPK.sbproj parameter FoldUptake = 0.01
    fold_efflux               <- fixed(0.01); label("Active efflux transporter clearance (iw to ew), as a fold of tissue plasma flow")  # Section 4.5 / supplement BasePBPK.sbproj parameter FoldEfflux = 0.01

    # Metabolic clearance is resolved PER COMPARTMENT because the paper
    # varies it both globally ("metabolism in all tissues", Tables 1, 3 and
    # 4) and in the liver alone ("metabolism in liver", Table 2), and
    # because Section 4.1 states that "metabolism can occur in any
    # subcompartment". Each is a fold-multiple of that compartment's
    # plasma flow, except the blood-cell clearances, which are a fold of
    # its blood-cell flow (supplement BasePBPK.sbproj rules, e.g.
    # AdiposeCLintiw = FoldClearanceIW * Qplasmaart2adipose and
    # AdiposeCLintrbc = FoldClearanceRBC * Qrbcart2adipose). Setting all
    # twelve tissues to one value reproduces the supplement's single
    # global multiplier exactly.
    fold_clint_iw_adipose     <- fixed(0.01); label("Metabolic clearance in intracellular water of adipose, as a fold of its plasma flow")  # Section 4.4.2 / supplement BasePBPK.sbproj parameter FoldClearanceIW = 0.01
    fold_clint_iw_bone        <- fixed(0.01); label("Metabolic clearance in intracellular water of bone, as a fold of its plasma flow")  # Section 4.4.2 / supplement BasePBPK.sbproj parameter FoldClearanceIW = 0.01
    fold_clint_iw_brain       <- fixed(0.01); label("Metabolic clearance in intracellular water of brain, as a fold of its plasma flow")  # Section 4.4.2 / supplement BasePBPK.sbproj parameter FoldClearanceIW = 0.01
    fold_clint_iw_heart       <- fixed(0.01); label("Metabolic clearance in intracellular water of heart, as a fold of its plasma flow")  # Section 4.4.2 / supplement BasePBPK.sbproj parameter FoldClearanceIW = 0.01
    fold_clint_iw_kidney      <- fixed(0.01); label("Metabolic clearance in intracellular water of kidney, as a fold of its plasma flow")  # Section 4.4.2 / supplement BasePBPK.sbproj parameter FoldClearanceIW = 0.01
    fold_clint_iw_muscle      <- fixed(0.01); label("Metabolic clearance in intracellular water of muscle, as a fold of its plasma flow")  # Section 4.4.2 / supplement BasePBPK.sbproj parameter FoldClearanceIW = 0.01
    fold_clint_iw_skin        <- fixed(0.01); label("Metabolic clearance in intracellular water of skin, as a fold of its plasma flow")  # Section 4.4.2 / supplement BasePBPK.sbproj parameter FoldClearanceIW = 0.01
    fold_clint_iw_liver       <- fixed(0.01); label("Metabolic clearance in intracellular water of liver, as a fold of its plasma flow")  # Section 4.4.2 / supplement BasePBPK.sbproj parameter FoldClearanceIW = 0.01
    fold_clint_iw_pancreas    <- fixed(0.01); label("Metabolic clearance in intracellular water of pancreas, as a fold of its plasma flow")  # Section 4.4.2 / supplement BasePBPK.sbproj parameter FoldClearanceIW = 0.01
    fold_clint_iw_spleen      <- fixed(0.01); label("Metabolic clearance in intracellular water of spleen, as a fold of its plasma flow")  # Section 4.4.2 / supplement BasePBPK.sbproj parameter FoldClearanceIW = 0.01
    fold_clint_iw_gut         <- fixed(0.01); label("Metabolic clearance in intracellular water of gut, as a fold of its plasma flow")  # Section 4.4.2 / supplement BasePBPK.sbproj parameter FoldClearanceIW = 0.01
    fold_clint_iw_lung        <- fixed(0.01); label("Metabolic clearance in intracellular water of lung, as a fold of its plasma flow")  # Section 4.4.2 / supplement BasePBPK.sbproj parameter FoldClearanceIW = 0.01
    fold_clint_ew_adipose     <- fixed(0.01); label("Metabolic clearance in extracellular water of adipose, as a fold of its plasma flow")  # Section 4.4.2 / supplement BasePBPK.sbproj parameter FoldClearanceEW = 0.01
    fold_clint_ew_bone        <- fixed(0.01); label("Metabolic clearance in extracellular water of bone, as a fold of its plasma flow")  # Section 4.4.2 / supplement BasePBPK.sbproj parameter FoldClearanceEW = 0.01
    fold_clint_ew_brain       <- fixed(0.01); label("Metabolic clearance in extracellular water of brain, as a fold of its plasma flow")  # Section 4.4.2 / supplement BasePBPK.sbproj parameter FoldClearanceEW = 0.01
    fold_clint_ew_heart       <- fixed(0.01); label("Metabolic clearance in extracellular water of heart, as a fold of its plasma flow")  # Section 4.4.2 / supplement BasePBPK.sbproj parameter FoldClearanceEW = 0.01
    fold_clint_ew_kidney      <- fixed(0.01); label("Metabolic clearance in extracellular water of kidney, as a fold of its plasma flow")  # Section 4.4.2 / supplement BasePBPK.sbproj parameter FoldClearanceEW = 0.01
    fold_clint_ew_muscle      <- fixed(0.01); label("Metabolic clearance in extracellular water of muscle, as a fold of its plasma flow")  # Section 4.4.2 / supplement BasePBPK.sbproj parameter FoldClearanceEW = 0.01
    fold_clint_ew_skin        <- fixed(0.01); label("Metabolic clearance in extracellular water of skin, as a fold of its plasma flow")  # Section 4.4.2 / supplement BasePBPK.sbproj parameter FoldClearanceEW = 0.01
    fold_clint_ew_liver       <- fixed(0.01); label("Metabolic clearance in extracellular water of liver, as a fold of its plasma flow")  # Section 4.4.2 / supplement BasePBPK.sbproj parameter FoldClearanceEW = 0.01
    fold_clint_ew_pancreas    <- fixed(0.01); label("Metabolic clearance in extracellular water of pancreas, as a fold of its plasma flow")  # Section 4.4.2 / supplement BasePBPK.sbproj parameter FoldClearanceEW = 0.01
    fold_clint_ew_spleen      <- fixed(0.01); label("Metabolic clearance in extracellular water of spleen, as a fold of its plasma flow")  # Section 4.4.2 / supplement BasePBPK.sbproj parameter FoldClearanceEW = 0.01
    fold_clint_ew_gut         <- fixed(0.01); label("Metabolic clearance in extracellular water of gut, as a fold of its plasma flow")  # Section 4.4.2 / supplement BasePBPK.sbproj parameter FoldClearanceEW = 0.01
    fold_clint_ew_lung        <- fixed(0.01); label("Metabolic clearance in extracellular water of lung, as a fold of its plasma flow")  # Section 4.4.2 / supplement BasePBPK.sbproj parameter FoldClearanceEW = 0.01
    fold_clint_plasma_adipose <- fixed(0.01); label("Metabolic clearance in residual plasma of adipose, as a fold of its plasma flow")  # Section 4.4.2 / supplement BasePBPK.sbproj parameter FoldClearancePLASMA = 0.01
    fold_clint_plasma_bone    <- fixed(0.01); label("Metabolic clearance in residual plasma of bone, as a fold of its plasma flow")  # Section 4.4.2 / supplement BasePBPK.sbproj parameter FoldClearancePLASMA = 0.01
    fold_clint_plasma_brain   <- fixed(0.01); label("Metabolic clearance in residual plasma of brain, as a fold of its plasma flow")  # Section 4.4.2 / supplement BasePBPK.sbproj parameter FoldClearancePLASMA = 0.01
    fold_clint_plasma_heart   <- fixed(0.01); label("Metabolic clearance in residual plasma of heart, as a fold of its plasma flow")  # Section 4.4.2 / supplement BasePBPK.sbproj parameter FoldClearancePLASMA = 0.01
    fold_clint_plasma_kidney  <- fixed(0.01); label("Metabolic clearance in residual plasma of kidney, as a fold of its plasma flow")  # Section 4.4.2 / supplement BasePBPK.sbproj parameter FoldClearancePLASMA = 0.01
    fold_clint_plasma_muscle  <- fixed(0.01); label("Metabolic clearance in residual plasma of muscle, as a fold of its plasma flow")  # Section 4.4.2 / supplement BasePBPK.sbproj parameter FoldClearancePLASMA = 0.01
    fold_clint_plasma_skin    <- fixed(0.01); label("Metabolic clearance in residual plasma of skin, as a fold of its plasma flow")  # Section 4.4.2 / supplement BasePBPK.sbproj parameter FoldClearancePLASMA = 0.01
    fold_clint_plasma_liver   <- fixed(0.01); label("Metabolic clearance in residual plasma of liver, as a fold of its plasma flow")  # Section 4.4.2 / supplement BasePBPK.sbproj parameter FoldClearancePLASMA = 0.01
    fold_clint_plasma_pancreas <- fixed(0.01); label("Metabolic clearance in residual plasma of pancreas, as a fold of its plasma flow")  # Section 4.4.2 / supplement BasePBPK.sbproj parameter FoldClearancePLASMA = 0.01
    fold_clint_plasma_spleen  <- fixed(0.01); label("Metabolic clearance in residual plasma of spleen, as a fold of its plasma flow")  # Section 4.4.2 / supplement BasePBPK.sbproj parameter FoldClearancePLASMA = 0.01
    fold_clint_plasma_gut     <- fixed(0.01); label("Metabolic clearance in residual plasma of gut, as a fold of its plasma flow")  # Section 4.4.2 / supplement BasePBPK.sbproj parameter FoldClearancePLASMA = 0.01
    fold_clint_plasma_lung    <- fixed(0.01); label("Metabolic clearance in residual plasma of lung, as a fold of its plasma flow")  # Section 4.4.2 / supplement BasePBPK.sbproj parameter FoldClearancePLASMA = 0.01
    fold_clint_plasma_venous  <- fixed(0.01); label("Metabolic clearance in residual plasma of venous, as a fold of its plasma flow")  # Section 4.4.2 / supplement BasePBPK.sbproj parameter FoldClearancePLASMA = 0.01
    fold_clint_plasma_arterial <- fixed(0.01); label("Metabolic clearance in residual plasma of arterial, as a fold of its plasma flow")  # Section 4.4.2 / supplement BasePBPK.sbproj parameter FoldClearancePLASMA = 0.01
    fold_clint_plasma_portal  <- fixed(0.01); label("Metabolic clearance in residual plasma of portal, as a fold of its plasma flow")  # Section 4.4.2 / supplement BasePBPK.sbproj parameter FoldClearancePLASMA = 0.01
    fold_clint_bc_adipose     <- fixed(0.01); label("Metabolic clearance in residual blood cells of adipose, as a fold of its blood-cell flow")  # Section 4.4.2 / supplement BasePBPK.sbproj parameter FoldClearanceRBC = 0.01
    fold_clint_bc_bone        <- fixed(0.01); label("Metabolic clearance in residual blood cells of bone, as a fold of its blood-cell flow")  # Section 4.4.2 / supplement BasePBPK.sbproj parameter FoldClearanceRBC = 0.01
    fold_clint_bc_brain       <- fixed(0.01); label("Metabolic clearance in residual blood cells of brain, as a fold of its blood-cell flow")  # Section 4.4.2 / supplement BasePBPK.sbproj parameter FoldClearanceRBC = 0.01
    fold_clint_bc_heart       <- fixed(0.01); label("Metabolic clearance in residual blood cells of heart, as a fold of its blood-cell flow")  # Section 4.4.2 / supplement BasePBPK.sbproj parameter FoldClearanceRBC = 0.01
    fold_clint_bc_kidney      <- fixed(0.01); label("Metabolic clearance in residual blood cells of kidney, as a fold of its blood-cell flow")  # Section 4.4.2 / supplement BasePBPK.sbproj parameter FoldClearanceRBC = 0.01
    fold_clint_bc_muscle      <- fixed(0.01); label("Metabolic clearance in residual blood cells of muscle, as a fold of its blood-cell flow")  # Section 4.4.2 / supplement BasePBPK.sbproj parameter FoldClearanceRBC = 0.01
    fold_clint_bc_skin        <- fixed(0.01); label("Metabolic clearance in residual blood cells of skin, as a fold of its blood-cell flow")  # Section 4.4.2 / supplement BasePBPK.sbproj parameter FoldClearanceRBC = 0.01
    fold_clint_bc_liver       <- fixed(0.01); label("Metabolic clearance in residual blood cells of liver, as a fold of its blood-cell flow")  # Section 4.4.2 / supplement BasePBPK.sbproj parameter FoldClearanceRBC = 0.01
    fold_clint_bc_pancreas    <- fixed(0.01); label("Metabolic clearance in residual blood cells of pancreas, as a fold of its blood-cell flow")  # Section 4.4.2 / supplement BasePBPK.sbproj parameter FoldClearanceRBC = 0.01
    fold_clint_bc_spleen      <- fixed(0.01); label("Metabolic clearance in residual blood cells of spleen, as a fold of its blood-cell flow")  # Section 4.4.2 / supplement BasePBPK.sbproj parameter FoldClearanceRBC = 0.01
    fold_clint_bc_gut         <- fixed(0.01); label("Metabolic clearance in residual blood cells of gut, as a fold of its blood-cell flow")  # Section 4.4.2 / supplement BasePBPK.sbproj parameter FoldClearanceRBC = 0.01
    fold_clint_bc_lung        <- fixed(0.01); label("Metabolic clearance in residual blood cells of lung, as a fold of its blood-cell flow")  # Section 4.4.2 / supplement BasePBPK.sbproj parameter FoldClearanceRBC = 0.01
    fold_clint_bc_venous      <- fixed(0.01); label("Metabolic clearance in residual blood cells of venous, as a fold of its blood-cell flow")  # Section 4.4.2 / supplement BasePBPK.sbproj parameter FoldClearanceRBC = 0.01
    fold_clint_bc_arterial    <- fixed(0.01); label("Metabolic clearance in residual blood cells of arterial, as a fold of its blood-cell flow")  # Section 4.4.2 / supplement BasePBPK.sbproj parameter FoldClearanceRBC = 0.01
    fold_clint_bc_portal      <- fixed(0.01); label("Metabolic clearance in residual blood cells of portal, as a fold of its blood-cell flow")  # Section 4.4.2 / supplement BasePBPK.sbproj parameter FoldClearanceRBC = 0.01

    # Paper Section 4.3.2: no binding and no ionisation is assumed, so
    # fu_tc = fu_tp = fu_ew = fu_iw = 1 and fi_tc = fi_tp = fi_ew = fi_iw
    # = 1 in all tissues and blood.
    fu_bc     <- fixed(1); label("Unbound fraction in blood cells")            # Section 4.3.2 (fu_tc = 1)
    fu_plasma <- fixed(1); label("Unbound fraction in plasma")                 # Section 4.3.2 (fu_tp = 1)
    fu_ew     <- fixed(1); label("Unbound fraction in extracellular water")    # Section 4.3.2 (fu_ew = 1)
    fu_iw     <- fixed(1); label("Unbound fraction in intracellular water")    # Section 4.3.2 (fu_iw = 1)
    fi_bc     <- fixed(1); label("Unionised fraction in blood cells")          # Section 4.3.2 (fi_tc = 1)
    fi_plasma <- fixed(1); label("Unionised fraction in plasma")               # Section 4.3.2 (fi_tp = 1)
    fi_ew     <- fixed(1); label("Unionised fraction in extracellular water")  # Section 4.3.2 (fi_ew = 1)
    fi_iw     <- fixed(1); label("Unionised fraction in intracellular water")  # Section 4.3.2 (fi_iw = 1)

    # Paper Section 4.3.3: ka = 1 1/h for fast oral absorption, 0.01 1/h
    # for slow oral absorption. The fast value is shipped as the default.
    lka <- fixed(log(1)); label("Oral absorption rate constant (1/h)")  # Section 4.3.3 (Ka = 1 1/h for fast oral absorption or Ka = 0.01 1/h for slow)

    # Numerical guard on the venous-plasma concentration in the Kp and Vd
    # denominators, so that both are finite before the first dose.
    cguard <- fixed(1e-60); label("Guard concentration added to the Kp / Vd denominator (mg/L)")  # supplement BasePBPK.sbproj parameter minorconcentration = 1e-60
  })

  model({
    # ============== 0. Individual parameters ==============
    ka <- exp(lka)

    # ============== 1. Subcompartment volumes (L) ==============
    # Density is 1 kg/L, so a tissue's volume in L equals its share of
    # body weight in kg. Within a tissue the residual blood splits into
    # blood cells and plasma by haematocrit, extracellular water is a
    # fixed fraction of the tissue, and intracellular water is whatever
    # is left (supplement BasePBPK.sbproj volume rules).
    v_blood          <- fvol_blood * bw / density
    v_venous_plasma  <- fven_blood * (1 - hct) * v_blood
    v_venous_bc      <- fven_blood * hct * v_blood
    v_arterial_plasma <- (1 - fven_blood) * (1 - hct) * v_blood
    v_arterial_bc     <- (1 - fven_blood) * hct * v_blood
    v_portal_plasma   <- (1 - hct) * v_portal
    v_portal_bc       <- hct * v_portal

    v_adipose            <- fvol_adipose * bw / density
    v_adipose_bc         <- hct * frb_adipose * v_adipose
    v_adipose_plasma     <- (1 - hct) * frb_adipose * v_adipose
    v_adipose_ew         <- few_adipose * v_adipose
    v_adipose_iw         <- v_adipose - (v_adipose_bc + v_adipose_plasma + v_adipose_ew)
    v_bone               <- fvol_bone * bw / density
    v_bone_bc            <- hct * frb_bone * v_bone
    v_bone_plasma        <- (1 - hct) * frb_bone * v_bone
    v_bone_ew            <- few_bone * v_bone
    v_bone_iw            <- v_bone - (v_bone_bc + v_bone_plasma + v_bone_ew)
    v_brain              <- fvol_brain * bw / density
    v_brain_bc           <- hct * frb_brain * v_brain
    v_brain_plasma       <- (1 - hct) * frb_brain * v_brain
    v_brain_ew           <- few_brain * v_brain
    v_brain_iw           <- v_brain - (v_brain_bc + v_brain_plasma + v_brain_ew)
    v_heart              <- fvol_heart * bw / density
    v_heart_bc           <- hct * frb_heart * v_heart
    v_heart_plasma       <- (1 - hct) * frb_heart * v_heart
    v_heart_ew           <- few_heart * v_heart
    v_heart_iw           <- v_heart - (v_heart_bc + v_heart_plasma + v_heart_ew)
    v_kidney             <- fvol_kidney * bw / density
    v_kidney_bc          <- hct * frb_kidney * v_kidney
    v_kidney_plasma      <- (1 - hct) * frb_kidney * v_kidney
    v_kidney_ew          <- few_kidney * v_kidney
    v_kidney_iw          <- v_kidney - (v_kidney_bc + v_kidney_plasma + v_kidney_ew)
    v_muscle             <- fvol_muscle * bw / density
    v_muscle_bc          <- hct * frb_muscle * v_muscle
    v_muscle_plasma      <- (1 - hct) * frb_muscle * v_muscle
    v_muscle_ew          <- few_muscle * v_muscle
    v_muscle_iw          <- v_muscle - (v_muscle_bc + v_muscle_plasma + v_muscle_ew)
    v_skin               <- fvol_skin * bw / density
    v_skin_bc            <- hct * frb_skin * v_skin
    v_skin_plasma        <- (1 - hct) * frb_skin * v_skin
    v_skin_ew            <- few_skin * v_skin
    v_skin_iw            <- v_skin - (v_skin_bc + v_skin_plasma + v_skin_ew)
    v_liver              <- fvol_liver * bw / density
    v_liver_bc           <- hct * frb_liver * v_liver
    v_liver_plasma       <- (1 - hct) * frb_liver * v_liver
    v_liver_ew           <- few_liver * v_liver
    v_liver_iw           <- v_liver - (v_liver_bc + v_liver_plasma + v_liver_ew)
    v_pancreas           <- fvol_pancreas * bw / density
    v_pancreas_bc        <- hct * frb_pancreas * v_pancreas
    v_pancreas_plasma    <- (1 - hct) * frb_pancreas * v_pancreas
    v_pancreas_ew        <- few_pancreas * v_pancreas
    v_pancreas_iw        <- v_pancreas - (v_pancreas_bc + v_pancreas_plasma + v_pancreas_ew)
    v_spleen             <- fvol_spleen * bw / density
    v_spleen_bc          <- hct * frb_spleen * v_spleen
    v_spleen_plasma      <- (1 - hct) * frb_spleen * v_spleen
    v_spleen_ew          <- few_spleen * v_spleen
    v_spleen_iw          <- v_spleen - (v_spleen_bc + v_spleen_plasma + v_spleen_ew)
    v_gut                <- fvol_gut * bw / density
    v_gut_bc             <- hct * frb_gut * v_gut
    v_gut_plasma         <- (1 - hct) * frb_gut * v_gut
    v_gut_ew             <- few_gut * v_gut
    v_gut_iw             <- v_gut - (v_gut_bc + v_gut_plasma + v_gut_ew)
    v_lung               <- fvol_lung * bw / density
    v_lung_bc            <- hct * frb_lung * v_lung
    v_lung_plasma        <- (1 - hct) * frb_lung * v_lung
    v_lung_ew            <- few_lung * v_lung
    v_lung_iw            <- v_lung - (v_lung_bc + v_lung_plasma + v_lung_ew)

    # ============== 2. Plasma and blood-cell flows (L/h) ==============
    # Blood flow splits into a plasma flow and a blood-cell flow by
    # haematocrit (paper Section 4.1: "blood flow consists of plasma flow
    # (Qtp) and blood cell flow (Qtc)").
    qplasma_total <- (1 - hct) * qc
    qbc_total     <- hct * qc
    qplasma_adipose      <- fq_adipose * qplasma_total
    qbc_adipose          <- fq_adipose * qbc_total
    qplasma_bone         <- fq_bone * qplasma_total
    qbc_bone             <- fq_bone * qbc_total
    qplasma_brain        <- fq_brain * qplasma_total
    qbc_brain            <- fq_brain * qbc_total
    qplasma_heart        <- fq_heart * qplasma_total
    qbc_heart            <- fq_heart * qbc_total
    qplasma_kidney       <- fq_kidney * qplasma_total
    qbc_kidney           <- fq_kidney * qbc_total
    qplasma_muscle       <- fq_muscle * qplasma_total
    qbc_muscle           <- fq_muscle * qbc_total
    qplasma_skin         <- fq_skin * qplasma_total
    qbc_skin             <- fq_skin * qbc_total
    qplasma_pancreas     <- fq_pancreas * qplasma_total
    qbc_pancreas         <- fq_pancreas * qbc_total
    qplasma_spleen       <- fq_spleen * qplasma_total
    qbc_spleen           <- fq_spleen * qbc_total
    qplasma_gut          <- fq_gut * qplasma_total
    qbc_gut              <- fq_gut * qbc_total
    # Lung is perfused by the whole cardiac output (Table A1: 100%).
    qplasma_lung  <- qplasma_total
    qbc_lung      <- qbc_total
    # Liver: hepatic artery plus portal vein (paper Eqs 25-26, 29-30).
    qplasma_liver_art <- fq_liver_arterial * qplasma_total
    qbc_liver_art     <- fq_liver_arterial * qbc_total
    qplasma_portal    <- qplasma_pancreas + qplasma_spleen + qplasma_gut
    qbc_portal        <- qbc_pancreas + qbc_spleen + qbc_gut
    qplasma_liver     <- qplasma_liver_art + qplasma_portal
    qbc_liver         <- qbc_liver_art + qbc_portal
    # Arterial-to-venous shunt closing the Eq 33-34 flow balance. Under
    # Table A1 the eight venous-draining flows sum to exactly 100% of
    # cardiac output, so this term is exactly zero; it is retained so that
    # mass balance still holds if the flow fractions are changed.
    qplasma_shunt <- qplasma_total - (qplasma_adipose + qplasma_bone +
      qplasma_brain + qplasma_heart + qplasma_kidney + qplasma_muscle +
      qplasma_skin + qplasma_liver)
    qbc_shunt     <- qbc_total - (qbc_adipose + qbc_bone + qbc_brain +
      qbc_heart + qbc_kidney + qbc_muscle + qbc_skin + qbc_liver)

    # ============== 3. Subcompartment concentrations (mg/L) ==============
    # States hold amounts (mg); concentrations are amount / volume.
    c_venous_plasma   <- venous_plasma / v_venous_plasma
    c_venous_bc       <- venous_bc / v_venous_bc
    c_arterial_plasma <- arterial_plasma / v_arterial_plasma
    c_arterial_bc     <- arterial_bc / v_arterial_bc
    c_portal_plasma   <- portal_plasma / v_portal_plasma
    c_portal_bc       <- portal_bc / v_portal_bc
    c_adipose_bc         <- adipose_bc / v_adipose_bc
    c_adipose_plasma     <- adipose_plasma / v_adipose_plasma
    c_adipose_ew         <- adipose_ew / v_adipose_ew
    c_adipose_iw         <- adipose_iw / v_adipose_iw
    c_bone_bc            <- bone_bc / v_bone_bc
    c_bone_plasma        <- bone_plasma / v_bone_plasma
    c_bone_ew            <- bone_ew / v_bone_ew
    c_bone_iw            <- bone_iw / v_bone_iw
    c_brain_bc           <- brain_bc / v_brain_bc
    c_brain_plasma       <- brain_plasma / v_brain_plasma
    c_brain_ew           <- brain_ew / v_brain_ew
    c_brain_iw           <- brain_iw / v_brain_iw
    c_heart_bc           <- heart_bc / v_heart_bc
    c_heart_plasma       <- heart_plasma / v_heart_plasma
    c_heart_ew           <- heart_ew / v_heart_ew
    c_heart_iw           <- heart_iw / v_heart_iw
    c_kidney_bc          <- kidney_bc / v_kidney_bc
    c_kidney_plasma      <- kidney_plasma / v_kidney_plasma
    c_kidney_ew          <- kidney_ew / v_kidney_ew
    c_kidney_iw          <- kidney_iw / v_kidney_iw
    c_muscle_bc          <- muscle_bc / v_muscle_bc
    c_muscle_plasma      <- muscle_plasma / v_muscle_plasma
    c_muscle_ew          <- muscle_ew / v_muscle_ew
    c_muscle_iw          <- muscle_iw / v_muscle_iw
    c_skin_bc            <- skin_bc / v_skin_bc
    c_skin_plasma        <- skin_plasma / v_skin_plasma
    c_skin_ew            <- skin_ew / v_skin_ew
    c_skin_iw            <- skin_iw / v_skin_iw
    c_liver_bc           <- liver_bc / v_liver_bc
    c_liver_plasma       <- liver_plasma / v_liver_plasma
    c_liver_ew           <- liver_ew / v_liver_ew
    c_liver_iw           <- liver_iw / v_liver_iw
    c_pancreas_bc        <- pancreas_bc / v_pancreas_bc
    c_pancreas_plasma    <- pancreas_plasma / v_pancreas_plasma
    c_pancreas_ew        <- pancreas_ew / v_pancreas_ew
    c_pancreas_iw        <- pancreas_iw / v_pancreas_iw
    c_spleen_bc          <- spleen_bc / v_spleen_bc
    c_spleen_plasma      <- spleen_plasma / v_spleen_plasma
    c_spleen_ew          <- spleen_ew / v_spleen_ew
    c_spleen_iw          <- spleen_iw / v_spleen_iw
    c_gut_bc             <- gut_bc / v_gut_bc
    c_gut_plasma         <- gut_plasma / v_gut_plasma
    c_gut_ew             <- gut_ew / v_gut_ew
    c_gut_iw             <- gut_iw / v_gut_iw
    c_lung_bc            <- lung_bc / v_lung_bc
    c_lung_plasma        <- lung_plasma / v_lung_plasma
    c_lung_ew            <- lung_ew / v_lung_ew
    c_lung_iw            <- lung_iw / v_lung_iw

    # ============== 4. Permeability, clearance and flux terms ==============
    # PS and every clearance is a fold-multiple of that tissue's plasma
    # flow (blood-cell flow for the blood-cell metabolic clearance), per
    # paper Sections 4.4.1-4.4.3 and the supplement BasePBPK.sbproj rules
    # (e.g. AdiposePSew2iw = FoldEwIw * Qplasmaart2adipose). Fluxes follow
    # paper Eqs 8-16.

    # ---- adipose ----
    ps_adipose_plasma_ew <- fold_ps_plasma_ew * qplasma_adipose
    ps_adipose_ew_iw     <- fold_ps_ew_iw * qplasma_adipose
    cl_adipose_iw     <- fold_clint_iw_adipose * qplasma_adipose
    cl_adipose_ew     <- fold_clint_ew_adipose * qplasma_adipose
    cl_adipose_plasma <- fold_clint_plasma_adipose * qplasma_adipose
    cl_adipose_bc     <- fold_clint_bc_adipose * qbc_adipose
    cl_adipose_uptake <- fold_uptake * qplasma_adipose
    cl_adipose_efflux <- fold_efflux * qplasma_adipose
    j_adipose_bc_plasma <- ps_bc_plasma * (c_adipose_bc * fu_bc * fi_bc - c_adipose_plasma * fu_plasma * fi_plasma)
    j_adipose_plasma_ew <- ps_adipose_plasma_ew * (c_adipose_plasma * fu_plasma * fi_plasma - c_adipose_ew * fu_ew * fi_ew)
    j_adipose_ew_iw     <- ps_adipose_ew_iw * (c_adipose_ew * fu_ew * fi_ew - c_adipose_iw * fu_iw * fi_iw)
    j_adipose_uptake    <- cl_adipose_uptake * c_adipose_ew * fu_ew
    j_adipose_efflux    <- cl_adipose_efflux * c_adipose_iw * fu_iw
    j_adipose_met_iw     <- cl_adipose_iw * c_adipose_iw * fu_iw
    j_adipose_met_ew     <- cl_adipose_ew * c_adipose_ew * fu_ew
    j_adipose_met_plasma <- cl_adipose_plasma * c_adipose_plasma * fu_plasma
    j_adipose_met_bc     <- cl_adipose_bc * c_adipose_bc * fu_bc

    # ---- bone ----
    ps_bone_plasma_ew <- fold_ps_plasma_ew * qplasma_bone
    ps_bone_ew_iw     <- fold_ps_ew_iw * qplasma_bone
    cl_bone_iw     <- fold_clint_iw_bone * qplasma_bone
    cl_bone_ew     <- fold_clint_ew_bone * qplasma_bone
    cl_bone_plasma <- fold_clint_plasma_bone * qplasma_bone
    cl_bone_bc     <- fold_clint_bc_bone * qbc_bone
    cl_bone_uptake <- fold_uptake * qplasma_bone
    cl_bone_efflux <- fold_efflux * qplasma_bone
    j_bone_bc_plasma <- ps_bc_plasma * (c_bone_bc * fu_bc * fi_bc - c_bone_plasma * fu_plasma * fi_plasma)
    j_bone_plasma_ew <- ps_bone_plasma_ew * (c_bone_plasma * fu_plasma * fi_plasma - c_bone_ew * fu_ew * fi_ew)
    j_bone_ew_iw     <- ps_bone_ew_iw * (c_bone_ew * fu_ew * fi_ew - c_bone_iw * fu_iw * fi_iw)
    j_bone_uptake    <- cl_bone_uptake * c_bone_ew * fu_ew
    j_bone_efflux    <- cl_bone_efflux * c_bone_iw * fu_iw
    j_bone_met_iw     <- cl_bone_iw * c_bone_iw * fu_iw
    j_bone_met_ew     <- cl_bone_ew * c_bone_ew * fu_ew
    j_bone_met_plasma <- cl_bone_plasma * c_bone_plasma * fu_plasma
    j_bone_met_bc     <- cl_bone_bc * c_bone_bc * fu_bc

    # ---- brain ----
    ps_brain_plasma_ew <- fold_ps_plasma_ew * qplasma_brain
    ps_brain_ew_iw     <- fold_ps_ew_iw * qplasma_brain
    cl_brain_iw     <- fold_clint_iw_brain * qplasma_brain
    cl_brain_ew     <- fold_clint_ew_brain * qplasma_brain
    cl_brain_plasma <- fold_clint_plasma_brain * qplasma_brain
    cl_brain_bc     <- fold_clint_bc_brain * qbc_brain
    cl_brain_uptake <- fold_uptake * qplasma_brain
    cl_brain_efflux <- fold_efflux * qplasma_brain
    j_brain_bc_plasma <- ps_bc_plasma * (c_brain_bc * fu_bc * fi_bc - c_brain_plasma * fu_plasma * fi_plasma)
    j_brain_plasma_ew <- ps_brain_plasma_ew * (c_brain_plasma * fu_plasma * fi_plasma - c_brain_ew * fu_ew * fi_ew)
    j_brain_ew_iw     <- ps_brain_ew_iw * (c_brain_ew * fu_ew * fi_ew - c_brain_iw * fu_iw * fi_iw)
    j_brain_uptake    <- cl_brain_uptake * c_brain_ew * fu_ew
    j_brain_efflux    <- cl_brain_efflux * c_brain_iw * fu_iw
    j_brain_met_iw     <- cl_brain_iw * c_brain_iw * fu_iw
    j_brain_met_ew     <- cl_brain_ew * c_brain_ew * fu_ew
    j_brain_met_plasma <- cl_brain_plasma * c_brain_plasma * fu_plasma
    j_brain_met_bc     <- cl_brain_bc * c_brain_bc * fu_bc

    # ---- heart ----
    ps_heart_plasma_ew <- fold_ps_plasma_ew * qplasma_heart
    ps_heart_ew_iw     <- fold_ps_ew_iw * qplasma_heart
    cl_heart_iw     <- fold_clint_iw_heart * qplasma_heart
    cl_heart_ew     <- fold_clint_ew_heart * qplasma_heart
    cl_heart_plasma <- fold_clint_plasma_heart * qplasma_heart
    cl_heart_bc     <- fold_clint_bc_heart * qbc_heart
    cl_heart_uptake <- fold_uptake * qplasma_heart
    cl_heart_efflux <- fold_efflux * qplasma_heart
    j_heart_bc_plasma <- ps_bc_plasma * (c_heart_bc * fu_bc * fi_bc - c_heart_plasma * fu_plasma * fi_plasma)
    j_heart_plasma_ew <- ps_heart_plasma_ew * (c_heart_plasma * fu_plasma * fi_plasma - c_heart_ew * fu_ew * fi_ew)
    j_heart_ew_iw     <- ps_heart_ew_iw * (c_heart_ew * fu_ew * fi_ew - c_heart_iw * fu_iw * fi_iw)
    j_heart_uptake    <- cl_heart_uptake * c_heart_ew * fu_ew
    j_heart_efflux    <- cl_heart_efflux * c_heart_iw * fu_iw
    j_heart_met_iw     <- cl_heart_iw * c_heart_iw * fu_iw
    j_heart_met_ew     <- cl_heart_ew * c_heart_ew * fu_ew
    j_heart_met_plasma <- cl_heart_plasma * c_heart_plasma * fu_plasma
    j_heart_met_bc     <- cl_heart_bc * c_heart_bc * fu_bc

    # ---- kidney ----
    ps_kidney_plasma_ew <- fold_ps_plasma_ew * qplasma_kidney
    ps_kidney_ew_iw     <- fold_ps_ew_iw * qplasma_kidney
    cl_kidney_iw     <- fold_clint_iw_kidney * qplasma_kidney
    cl_kidney_ew     <- fold_clint_ew_kidney * qplasma_kidney
    cl_kidney_plasma <- fold_clint_plasma_kidney * qplasma_kidney
    cl_kidney_bc     <- fold_clint_bc_kidney * qbc_kidney
    cl_kidney_uptake <- fold_uptake * qplasma_kidney
    cl_kidney_efflux <- fold_efflux * qplasma_kidney
    j_kidney_bc_plasma <- ps_bc_plasma * (c_kidney_bc * fu_bc * fi_bc - c_kidney_plasma * fu_plasma * fi_plasma)
    j_kidney_plasma_ew <- ps_kidney_plasma_ew * (c_kidney_plasma * fu_plasma * fi_plasma - c_kidney_ew * fu_ew * fi_ew)
    j_kidney_ew_iw     <- ps_kidney_ew_iw * (c_kidney_ew * fu_ew * fi_ew - c_kidney_iw * fu_iw * fi_iw)
    j_kidney_uptake    <- cl_kidney_uptake * c_kidney_ew * fu_ew
    j_kidney_efflux    <- cl_kidney_efflux * c_kidney_iw * fu_iw
    j_kidney_met_iw     <- cl_kidney_iw * c_kidney_iw * fu_iw
    j_kidney_met_ew     <- cl_kidney_ew * c_kidney_ew * fu_ew
    j_kidney_met_plasma <- cl_kidney_plasma * c_kidney_plasma * fu_plasma
    j_kidney_met_bc     <- cl_kidney_bc * c_kidney_bc * fu_bc

    # ---- muscle ----
    ps_muscle_plasma_ew <- fold_ps_plasma_ew * qplasma_muscle
    ps_muscle_ew_iw     <- fold_ps_ew_iw * qplasma_muscle
    cl_muscle_iw     <- fold_clint_iw_muscle * qplasma_muscle
    cl_muscle_ew     <- fold_clint_ew_muscle * qplasma_muscle
    cl_muscle_plasma <- fold_clint_plasma_muscle * qplasma_muscle
    cl_muscle_bc     <- fold_clint_bc_muscle * qbc_muscle
    cl_muscle_uptake <- fold_uptake * qplasma_muscle
    cl_muscle_efflux <- fold_efflux * qplasma_muscle
    j_muscle_bc_plasma <- ps_bc_plasma * (c_muscle_bc * fu_bc * fi_bc - c_muscle_plasma * fu_plasma * fi_plasma)
    j_muscle_plasma_ew <- ps_muscle_plasma_ew * (c_muscle_plasma * fu_plasma * fi_plasma - c_muscle_ew * fu_ew * fi_ew)
    j_muscle_ew_iw     <- ps_muscle_ew_iw * (c_muscle_ew * fu_ew * fi_ew - c_muscle_iw * fu_iw * fi_iw)
    j_muscle_uptake    <- cl_muscle_uptake * c_muscle_ew * fu_ew
    j_muscle_efflux    <- cl_muscle_efflux * c_muscle_iw * fu_iw
    j_muscle_met_iw     <- cl_muscle_iw * c_muscle_iw * fu_iw
    j_muscle_met_ew     <- cl_muscle_ew * c_muscle_ew * fu_ew
    j_muscle_met_plasma <- cl_muscle_plasma * c_muscle_plasma * fu_plasma
    j_muscle_met_bc     <- cl_muscle_bc * c_muscle_bc * fu_bc

    # ---- skin ----
    ps_skin_plasma_ew <- fold_ps_plasma_ew * qplasma_skin
    ps_skin_ew_iw     <- fold_ps_ew_iw * qplasma_skin
    cl_skin_iw     <- fold_clint_iw_skin * qplasma_skin
    cl_skin_ew     <- fold_clint_ew_skin * qplasma_skin
    cl_skin_plasma <- fold_clint_plasma_skin * qplasma_skin
    cl_skin_bc     <- fold_clint_bc_skin * qbc_skin
    cl_skin_uptake <- fold_uptake * qplasma_skin
    cl_skin_efflux <- fold_efflux * qplasma_skin
    j_skin_bc_plasma <- ps_bc_plasma * (c_skin_bc * fu_bc * fi_bc - c_skin_plasma * fu_plasma * fi_plasma)
    j_skin_plasma_ew <- ps_skin_plasma_ew * (c_skin_plasma * fu_plasma * fi_plasma - c_skin_ew * fu_ew * fi_ew)
    j_skin_ew_iw     <- ps_skin_ew_iw * (c_skin_ew * fu_ew * fi_ew - c_skin_iw * fu_iw * fi_iw)
    j_skin_uptake    <- cl_skin_uptake * c_skin_ew * fu_ew
    j_skin_efflux    <- cl_skin_efflux * c_skin_iw * fu_iw
    j_skin_met_iw     <- cl_skin_iw * c_skin_iw * fu_iw
    j_skin_met_ew     <- cl_skin_ew * c_skin_ew * fu_ew
    j_skin_met_plasma <- cl_skin_plasma * c_skin_plasma * fu_plasma
    j_skin_met_bc     <- cl_skin_bc * c_skin_bc * fu_bc

    # ---- liver ----
    ps_liver_plasma_ew <- fold_ps_plasma_ew * qplasma_liver
    ps_liver_ew_iw     <- fold_ps_ew_iw * qplasma_liver
    cl_liver_iw     <- fold_clint_iw_liver * qplasma_liver
    cl_liver_ew     <- fold_clint_ew_liver * qplasma_liver
    cl_liver_plasma <- fold_clint_plasma_liver * qplasma_liver
    cl_liver_bc     <- fold_clint_bc_liver * qbc_liver
    cl_liver_uptake <- fold_uptake * qplasma_liver
    cl_liver_efflux <- fold_efflux * qplasma_liver
    j_liver_bc_plasma <- ps_bc_plasma * (c_liver_bc * fu_bc * fi_bc - c_liver_plasma * fu_plasma * fi_plasma)
    j_liver_plasma_ew <- ps_liver_plasma_ew * (c_liver_plasma * fu_plasma * fi_plasma - c_liver_ew * fu_ew * fi_ew)
    j_liver_ew_iw     <- ps_liver_ew_iw * (c_liver_ew * fu_ew * fi_ew - c_liver_iw * fu_iw * fi_iw)
    j_liver_uptake    <- cl_liver_uptake * c_liver_ew * fu_ew
    j_liver_efflux    <- cl_liver_efflux * c_liver_iw * fu_iw
    j_liver_met_iw     <- cl_liver_iw * c_liver_iw * fu_iw
    j_liver_met_ew     <- cl_liver_ew * c_liver_ew * fu_ew
    j_liver_met_plasma <- cl_liver_plasma * c_liver_plasma * fu_plasma
    j_liver_met_bc     <- cl_liver_bc * c_liver_bc * fu_bc

    # ---- pancreas ----
    ps_pancreas_plasma_ew <- fold_ps_plasma_ew * qplasma_pancreas
    ps_pancreas_ew_iw     <- fold_ps_ew_iw * qplasma_pancreas
    cl_pancreas_iw     <- fold_clint_iw_pancreas * qplasma_pancreas
    cl_pancreas_ew     <- fold_clint_ew_pancreas * qplasma_pancreas
    cl_pancreas_plasma <- fold_clint_plasma_pancreas * qplasma_pancreas
    cl_pancreas_bc     <- fold_clint_bc_pancreas * qbc_pancreas
    cl_pancreas_uptake <- fold_uptake * qplasma_pancreas
    cl_pancreas_efflux <- fold_efflux * qplasma_pancreas
    j_pancreas_bc_plasma <- ps_bc_plasma * (c_pancreas_bc * fu_bc * fi_bc - c_pancreas_plasma * fu_plasma * fi_plasma)
    j_pancreas_plasma_ew <- ps_pancreas_plasma_ew * (c_pancreas_plasma * fu_plasma * fi_plasma - c_pancreas_ew * fu_ew * fi_ew)
    j_pancreas_ew_iw     <- ps_pancreas_ew_iw * (c_pancreas_ew * fu_ew * fi_ew - c_pancreas_iw * fu_iw * fi_iw)
    j_pancreas_uptake    <- cl_pancreas_uptake * c_pancreas_ew * fu_ew
    j_pancreas_efflux    <- cl_pancreas_efflux * c_pancreas_iw * fu_iw
    j_pancreas_met_iw     <- cl_pancreas_iw * c_pancreas_iw * fu_iw
    j_pancreas_met_ew     <- cl_pancreas_ew * c_pancreas_ew * fu_ew
    j_pancreas_met_plasma <- cl_pancreas_plasma * c_pancreas_plasma * fu_plasma
    j_pancreas_met_bc     <- cl_pancreas_bc * c_pancreas_bc * fu_bc

    # ---- spleen ----
    ps_spleen_plasma_ew <- fold_ps_plasma_ew * qplasma_spleen
    ps_spleen_ew_iw     <- fold_ps_ew_iw * qplasma_spleen
    cl_spleen_iw     <- fold_clint_iw_spleen * qplasma_spleen
    cl_spleen_ew     <- fold_clint_ew_spleen * qplasma_spleen
    cl_spleen_plasma <- fold_clint_plasma_spleen * qplasma_spleen
    cl_spleen_bc     <- fold_clint_bc_spleen * qbc_spleen
    cl_spleen_uptake <- fold_uptake * qplasma_spleen
    cl_spleen_efflux <- fold_efflux * qplasma_spleen
    j_spleen_bc_plasma <- ps_bc_plasma * (c_spleen_bc * fu_bc * fi_bc - c_spleen_plasma * fu_plasma * fi_plasma)
    j_spleen_plasma_ew <- ps_spleen_plasma_ew * (c_spleen_plasma * fu_plasma * fi_plasma - c_spleen_ew * fu_ew * fi_ew)
    j_spleen_ew_iw     <- ps_spleen_ew_iw * (c_spleen_ew * fu_ew * fi_ew - c_spleen_iw * fu_iw * fi_iw)
    j_spleen_uptake    <- cl_spleen_uptake * c_spleen_ew * fu_ew
    j_spleen_efflux    <- cl_spleen_efflux * c_spleen_iw * fu_iw
    j_spleen_met_iw     <- cl_spleen_iw * c_spleen_iw * fu_iw
    j_spleen_met_ew     <- cl_spleen_ew * c_spleen_ew * fu_ew
    j_spleen_met_plasma <- cl_spleen_plasma * c_spleen_plasma * fu_plasma
    j_spleen_met_bc     <- cl_spleen_bc * c_spleen_bc * fu_bc

    # ---- gut ----
    ps_gut_plasma_ew <- fold_ps_plasma_ew * qplasma_gut
    ps_gut_ew_iw     <- fold_ps_ew_iw * qplasma_gut
    cl_gut_iw     <- fold_clint_iw_gut * qplasma_gut
    cl_gut_ew     <- fold_clint_ew_gut * qplasma_gut
    cl_gut_plasma <- fold_clint_plasma_gut * qplasma_gut
    cl_gut_bc     <- fold_clint_bc_gut * qbc_gut
    cl_gut_uptake <- fold_uptake * qplasma_gut
    cl_gut_efflux <- fold_efflux * qplasma_gut
    j_gut_bc_plasma <- ps_bc_plasma * (c_gut_bc * fu_bc * fi_bc - c_gut_plasma * fu_plasma * fi_plasma)
    j_gut_plasma_ew <- ps_gut_plasma_ew * (c_gut_plasma * fu_plasma * fi_plasma - c_gut_ew * fu_ew * fi_ew)
    j_gut_ew_iw     <- ps_gut_ew_iw * (c_gut_ew * fu_ew * fi_ew - c_gut_iw * fu_iw * fi_iw)
    j_gut_uptake    <- cl_gut_uptake * c_gut_ew * fu_ew
    j_gut_efflux    <- cl_gut_efflux * c_gut_iw * fu_iw
    j_gut_met_iw     <- cl_gut_iw * c_gut_iw * fu_iw
    j_gut_met_ew     <- cl_gut_ew * c_gut_ew * fu_ew
    j_gut_met_plasma <- cl_gut_plasma * c_gut_plasma * fu_plasma
    j_gut_met_bc     <- cl_gut_bc * c_gut_bc * fu_bc

    # ---- lung ----
    ps_lung_plasma_ew <- fold_ps_plasma_ew * qplasma_lung
    ps_lung_ew_iw     <- fold_ps_ew_iw * qplasma_lung
    cl_lung_iw     <- fold_clint_iw_lung * qplasma_lung
    cl_lung_ew     <- fold_clint_ew_lung * qplasma_lung
    cl_lung_plasma <- fold_clint_plasma_lung * qplasma_lung
    cl_lung_bc     <- fold_clint_bc_lung * qbc_lung
    cl_lung_uptake <- fold_uptake * qplasma_lung
    cl_lung_efflux <- fold_efflux * qplasma_lung
    j_lung_bc_plasma <- ps_bc_plasma * (c_lung_bc * fu_bc * fi_bc - c_lung_plasma * fu_plasma * fi_plasma)
    j_lung_plasma_ew <- ps_lung_plasma_ew * (c_lung_plasma * fu_plasma * fi_plasma - c_lung_ew * fu_ew * fi_ew)
    j_lung_ew_iw     <- ps_lung_ew_iw * (c_lung_ew * fu_ew * fi_ew - c_lung_iw * fu_iw * fi_iw)
    j_lung_uptake    <- cl_lung_uptake * c_lung_ew * fu_ew
    j_lung_efflux    <- cl_lung_efflux * c_lung_iw * fu_iw
    j_lung_met_iw     <- cl_lung_iw * c_lung_iw * fu_iw
    j_lung_met_ew     <- cl_lung_ew * c_lung_ew * fu_ew
    j_lung_met_plasma <- cl_lung_plasma * c_lung_plasma * fu_plasma
    j_lung_met_bc     <- cl_lung_bc * c_lung_bc * fu_bc

    # ---- blood compartments: only the bc <-> plasma exchange and the
    # subcompartment metabolic clearances exist (paper Section 4.1).
    # Venous and arterial blood-cell clearances scale with the total
    # blood-cell flow, their plasma clearances with the total plasma flow;
    # the portal vein scales with the portal flow (paper Eqs 27-28, 31-32,
    # 35-36).
    cl_venous_plasma <- fold_clint_plasma_venous * qplasma_total
    cl_venous_bc     <- fold_clint_bc_venous * qbc_total
    j_venous_bc_plasma  <- ps_bc_plasma * (c_venous_bc * fu_bc * fi_bc - c_venous_plasma * fu_plasma * fi_plasma)
    j_venous_met_plasma <- cl_venous_plasma * c_venous_plasma * fu_plasma
    j_venous_met_bc     <- cl_venous_bc * c_venous_bc * fu_bc
    cl_arterial_plasma <- fold_clint_plasma_arterial * qplasma_total
    cl_arterial_bc     <- fold_clint_bc_arterial * qbc_total
    j_arterial_bc_plasma  <- ps_bc_plasma * (c_arterial_bc * fu_bc * fi_bc - c_arterial_plasma * fu_plasma * fi_plasma)
    j_arterial_met_plasma <- cl_arterial_plasma * c_arterial_plasma * fu_plasma
    j_arterial_met_bc     <- cl_arterial_bc * c_arterial_bc * fu_bc
    cl_portal_plasma <- fold_clint_plasma_portal * qplasma_portal
    cl_portal_bc     <- fold_clint_bc_portal * qbc_portal
    j_portal_bc_plasma  <- ps_bc_plasma * (c_portal_bc * fu_bc * fi_bc - c_portal_plasma * fu_plasma * fi_plasma)
    j_portal_met_plasma <- cl_portal_plasma * c_portal_plasma * fu_plasma
    j_portal_met_bc     <- cl_portal_bc * c_portal_bc * fu_bc

    # ============== 5. ODE system (55 states, all amounts in mg) ==============

    # ---- adipose ----
    # Eq 4: intracellular water
    d/dt(adipose_iw) <- j_adipose_ew_iw + j_adipose_uptake - j_adipose_efflux - j_adipose_met_iw
    # Eq 5: extracellular water
    d/dt(adipose_ew) <- j_adipose_plasma_ew - j_adipose_ew_iw - j_adipose_uptake +
      j_adipose_efflux - j_adipose_met_ew
    # Eq 6: residual plasma
    d/dt(adipose_plasma) <- qplasma_adipose * (c_arterial_plasma - c_adipose_plasma) + j_adipose_bc_plasma -
      j_adipose_plasma_ew - j_adipose_met_plasma
    # Eq 7: residual blood cells
    d/dt(adipose_bc) <- qbc_adipose * (c_arterial_bc - c_adipose_bc) - j_adipose_bc_plasma - j_adipose_met_bc

    # ---- bone ----
    # Eq 4: intracellular water
    d/dt(bone_iw) <- j_bone_ew_iw + j_bone_uptake - j_bone_efflux - j_bone_met_iw
    # Eq 5: extracellular water
    d/dt(bone_ew) <- j_bone_plasma_ew - j_bone_ew_iw - j_bone_uptake +
      j_bone_efflux - j_bone_met_ew
    # Eq 6: residual plasma
    d/dt(bone_plasma) <- qplasma_bone * (c_arterial_plasma - c_bone_plasma) + j_bone_bc_plasma -
      j_bone_plasma_ew - j_bone_met_plasma
    # Eq 7: residual blood cells
    d/dt(bone_bc) <- qbc_bone * (c_arterial_bc - c_bone_bc) - j_bone_bc_plasma - j_bone_met_bc

    # ---- brain ----
    # Eq 4: intracellular water
    d/dt(brain_iw) <- j_brain_ew_iw + j_brain_uptake - j_brain_efflux - j_brain_met_iw
    # Eq 5: extracellular water
    d/dt(brain_ew) <- j_brain_plasma_ew - j_brain_ew_iw - j_brain_uptake +
      j_brain_efflux - j_brain_met_ew
    # Eq 6: residual plasma
    d/dt(brain_plasma) <- qplasma_brain * (c_arterial_plasma - c_brain_plasma) + j_brain_bc_plasma -
      j_brain_plasma_ew - j_brain_met_plasma
    # Eq 7: residual blood cells
    d/dt(brain_bc) <- qbc_brain * (c_arterial_bc - c_brain_bc) - j_brain_bc_plasma - j_brain_met_bc

    # ---- heart ----
    # Eq 4: intracellular water
    d/dt(heart_iw) <- j_heart_ew_iw + j_heart_uptake - j_heart_efflux - j_heart_met_iw
    # Eq 5: extracellular water
    d/dt(heart_ew) <- j_heart_plasma_ew - j_heart_ew_iw - j_heart_uptake +
      j_heart_efflux - j_heart_met_ew
    # Eq 6: residual plasma
    d/dt(heart_plasma) <- qplasma_heart * (c_arterial_plasma - c_heart_plasma) + j_heart_bc_plasma -
      j_heart_plasma_ew - j_heart_met_plasma
    # Eq 7: residual blood cells
    d/dt(heart_bc) <- qbc_heart * (c_arterial_bc - c_heart_bc) - j_heart_bc_plasma - j_heart_met_bc

    # ---- kidney ----
    # Eq 4: intracellular water
    d/dt(kidney_iw) <- j_kidney_ew_iw + j_kidney_uptake - j_kidney_efflux - j_kidney_met_iw
    # Eq 5: extracellular water
    d/dt(kidney_ew) <- j_kidney_plasma_ew - j_kidney_ew_iw - j_kidney_uptake +
      j_kidney_efflux - j_kidney_met_ew
    # Eq 6: residual plasma
    d/dt(kidney_plasma) <- qplasma_kidney * (c_arterial_plasma - c_kidney_plasma) + j_kidney_bc_plasma -
      j_kidney_plasma_ew - j_kidney_met_plasma
    # Eq 7: residual blood cells
    d/dt(kidney_bc) <- qbc_kidney * (c_arterial_bc - c_kidney_bc) - j_kidney_bc_plasma - j_kidney_met_bc

    # ---- muscle ----
    # Eq 4: intracellular water
    d/dt(muscle_iw) <- j_muscle_ew_iw + j_muscle_uptake - j_muscle_efflux - j_muscle_met_iw
    # Eq 5: extracellular water
    d/dt(muscle_ew) <- j_muscle_plasma_ew - j_muscle_ew_iw - j_muscle_uptake +
      j_muscle_efflux - j_muscle_met_ew
    # Eq 6: residual plasma
    d/dt(muscle_plasma) <- qplasma_muscle * (c_arterial_plasma - c_muscle_plasma) + j_muscle_bc_plasma -
      j_muscle_plasma_ew - j_muscle_met_plasma
    # Eq 7: residual blood cells
    d/dt(muscle_bc) <- qbc_muscle * (c_arterial_bc - c_muscle_bc) - j_muscle_bc_plasma - j_muscle_met_bc

    # ---- skin ----
    # Eq 4: intracellular water
    d/dt(skin_iw) <- j_skin_ew_iw + j_skin_uptake - j_skin_efflux - j_skin_met_iw
    # Eq 5: extracellular water
    d/dt(skin_ew) <- j_skin_plasma_ew - j_skin_ew_iw - j_skin_uptake +
      j_skin_efflux - j_skin_met_ew
    # Eq 6: residual plasma
    d/dt(skin_plasma) <- qplasma_skin * (c_arterial_plasma - c_skin_plasma) + j_skin_bc_plasma -
      j_skin_plasma_ew - j_skin_met_plasma
    # Eq 7: residual blood cells
    d/dt(skin_bc) <- qbc_skin * (c_arterial_bc - c_skin_bc) - j_skin_bc_plasma - j_skin_met_bc

    # ---- pancreas ----
    # Eq 4: intracellular water
    d/dt(pancreas_iw) <- j_pancreas_ew_iw + j_pancreas_uptake - j_pancreas_efflux - j_pancreas_met_iw
    # Eq 5: extracellular water
    d/dt(pancreas_ew) <- j_pancreas_plasma_ew - j_pancreas_ew_iw - j_pancreas_uptake +
      j_pancreas_efflux - j_pancreas_met_ew
    # Eq 6: residual plasma
    d/dt(pancreas_plasma) <- qplasma_pancreas * (c_arterial_plasma - c_pancreas_plasma) + j_pancreas_bc_plasma -
      j_pancreas_plasma_ew - j_pancreas_met_plasma
    # Eq 7: residual blood cells
    d/dt(pancreas_bc) <- qbc_pancreas * (c_arterial_bc - c_pancreas_bc) - j_pancreas_bc_plasma - j_pancreas_met_bc

    # ---- spleen ----
    # Eq 4: intracellular water
    d/dt(spleen_iw) <- j_spleen_ew_iw + j_spleen_uptake - j_spleen_efflux - j_spleen_met_iw
    # Eq 5: extracellular water
    d/dt(spleen_ew) <- j_spleen_plasma_ew - j_spleen_ew_iw - j_spleen_uptake +
      j_spleen_efflux - j_spleen_met_ew
    # Eq 6: residual plasma
    d/dt(spleen_plasma) <- qplasma_spleen * (c_arterial_plasma - c_spleen_plasma) + j_spleen_bc_plasma -
      j_spleen_plasma_ew - j_spleen_met_plasma
    # Eq 7: residual blood cells
    d/dt(spleen_bc) <- qbc_spleen * (c_arterial_bc - c_spleen_bc) - j_spleen_bc_plasma - j_spleen_met_bc

    # ---- gut ----
    # Eq 4: intracellular water
    d/dt(gut_iw) <- j_gut_ew_iw + j_gut_uptake - j_gut_efflux - j_gut_met_iw +
      ka * depot
    # Eq 5: extracellular water
    d/dt(gut_ew) <- j_gut_plasma_ew - j_gut_ew_iw - j_gut_uptake +
      j_gut_efflux - j_gut_met_ew
    # Eq 6: residual plasma
    d/dt(gut_plasma) <- qplasma_gut * (c_arterial_plasma - c_gut_plasma) + j_gut_bc_plasma -
      j_gut_plasma_ew - j_gut_met_plasma
    # Eq 7: residual blood cells
    d/dt(gut_bc) <- qbc_gut * (c_arterial_bc - c_gut_bc) - j_gut_bc_plasma - j_gut_met_bc

    # ---- lung ----
    # Eq 17: intracellular water
    d/dt(lung_iw) <- j_lung_ew_iw + j_lung_uptake - j_lung_efflux - j_lung_met_iw
    # Eq 18: extracellular water
    d/dt(lung_ew) <- j_lung_plasma_ew - j_lung_ew_iw - j_lung_uptake +
      j_lung_efflux - j_lung_met_ew
    # Eq 19: residual plasma
    d/dt(lung_plasma) <- qplasma_lung * (c_venous_plasma - c_lung_plasma) + j_lung_bc_plasma -
      j_lung_plasma_ew - j_lung_met_plasma
    # Eq 20: residual blood cells
    d/dt(lung_bc) <- qbc_lung * (c_venous_bc - c_lung_bc) - j_lung_bc_plasma - j_lung_met_bc

    # ---- liver ----
    # Eq 21: intracellular water
    d/dt(liver_iw) <- j_liver_ew_iw + j_liver_uptake - j_liver_efflux - j_liver_met_iw
    # Eq 22: extracellular water
    d/dt(liver_ew) <- j_liver_plasma_ew - j_liver_ew_iw - j_liver_uptake +
      j_liver_efflux - j_liver_met_ew
    # Eq 23: residual plasma
    d/dt(liver_plasma) <- qplasma_liver_art * c_arterial_plasma +
      qplasma_portal * c_portal_plasma -
      qplasma_liver * c_liver_plasma + j_liver_bc_plasma -
      j_liver_plasma_ew - j_liver_met_plasma
    # Eq 24: residual blood cells
    d/dt(liver_bc) <- qbc_liver_art * c_arterial_bc +
      qbc_portal * c_portal_bc -
      qbc_liver * c_liver_bc - j_liver_bc_plasma - j_liver_met_bc

    # ---- portal vein (Eqs 27-28) ----
    d/dt(portal_plasma) <- qplasma_pancreas * c_pancreas_plasma +
      qplasma_spleen * c_spleen_plasma + qplasma_gut * c_gut_plasma -
      qplasma_portal * c_portal_plasma + j_portal_bc_plasma -
      j_portal_met_plasma
    d/dt(portal_bc) <- qbc_pancreas * c_pancreas_bc +
      qbc_spleen * c_spleen_bc + qbc_gut * c_gut_bc -
      qbc_portal * c_portal_bc - j_portal_bc_plasma - j_portal_met_bc

    # ---- venous blood (Eqs 31-32) ----
    # The sum runs over the tissues that drain directly to venous blood:
    # adipose, bone, brain, heart, kidney, muscle, skin and liver (paper
    # text after Eq 32). The shunt term is zero under Table A1.
    d/dt(venous_plasma) <- qplasma_adipose * c_adipose_plasma +
      qplasma_bone * c_bone_plasma + qplasma_brain * c_brain_plasma +
      qplasma_heart * c_heart_plasma + qplasma_kidney * c_kidney_plasma +
      qplasma_muscle * c_muscle_plasma + qplasma_skin * c_skin_plasma +
      qplasma_liver * c_liver_plasma +
      qplasma_shunt * c_arterial_plasma -
      qplasma_total * c_venous_plasma + j_venous_bc_plasma -
      j_venous_met_plasma
    d/dt(venous_bc) <- qbc_adipose * c_adipose_bc +
      qbc_bone * c_bone_bc + qbc_brain * c_brain_bc +
      qbc_heart * c_heart_bc + qbc_kidney * c_kidney_bc +
      qbc_muscle * c_muscle_bc + qbc_skin * c_skin_bc +
      qbc_liver * c_liver_bc + qbc_shunt * c_arterial_bc -
      qbc_total * c_venous_bc - j_venous_bc_plasma - j_venous_met_bc

    # ---- arterial blood (Eqs 35-36) ----
    d/dt(arterial_plasma) <- qplasma_lung * (c_lung_plasma - c_arterial_plasma) -
      qplasma_shunt * c_arterial_plasma + j_arterial_bc_plasma -
      j_arterial_met_plasma
    d/dt(arterial_bc) <- qbc_lung * (c_lung_bc - c_arterial_bc) -
      qbc_shunt * c_arterial_bc - j_arterial_bc_plasma - j_arterial_met_bc

    # ---- gut lumen depot (Eq 37) ----
    # Eq 37 states the oral input as Dose * Ka * exp(-Ka * t) into the gut
    # intracellular water. A first-order depot is the identical input
    # function and is how the supplement implements it (BasePBPK.sbproj
    # reaction Lumen.drug -> iw_gut.drug at rate Ka * Lumen.drug).
    d/dt(depot) <- -ka * depot

    # ============== 6. Tissue concentrations, Kp(t) and Vd(t) ==============
    # Eq 38: the whole-tissue concentration is the total amount in the four
    # subcompartments divided by the total tissue volume.
    c_adipose_tissue     <- (adipose_iw + adipose_ew + adipose_plasma + adipose_bc) / v_adipose
    c_bone_tissue        <- (bone_iw + bone_ew + bone_plasma + bone_bc) / v_bone
    c_brain_tissue       <- (brain_iw + brain_ew + brain_plasma + brain_bc) / v_brain
    c_heart_tissue       <- (heart_iw + heart_ew + heart_plasma + heart_bc) / v_heart
    c_kidney_tissue      <- (kidney_iw + kidney_ew + kidney_plasma + kidney_bc) / v_kidney
    c_muscle_tissue      <- (muscle_iw + muscle_ew + muscle_plasma + muscle_bc) / v_muscle
    c_skin_tissue        <- (skin_iw + skin_ew + skin_plasma + skin_bc) / v_skin
    c_liver_tissue       <- (liver_iw + liver_ew + liver_plasma + liver_bc) / v_liver
    c_pancreas_tissue    <- (pancreas_iw + pancreas_ew + pancreas_plasma + pancreas_bc) / v_pancreas
    c_spleen_tissue      <- (spleen_iw + spleen_ew + spleen_plasma + spleen_bc) / v_spleen
    c_gut_tissue         <- (gut_iw + gut_ew + gut_plasma + gut_bc) / v_gut
    c_lung_tissue        <- (lung_iw + lung_ew + lung_plasma + lung_bc) / v_lung

    # Eq 39: tissue/plasma partition coefficient, using the venous plasma
    # concentration as the surrogate for systemic plasma.
    Kp_adipose     <- c_adipose_tissue / (c_venous_plasma + cguard)
    Kp_bone        <- c_bone_tissue / (c_venous_plasma + cguard)
    Kp_brain       <- c_brain_tissue / (c_venous_plasma + cguard)
    Kp_heart       <- c_heart_tissue / (c_venous_plasma + cguard)
    Kp_kidney      <- c_kidney_tissue / (c_venous_plasma + cguard)
    Kp_muscle      <- c_muscle_tissue / (c_venous_plasma + cguard)
    Kp_skin        <- c_skin_tissue / (c_venous_plasma + cguard)
    Kp_liver       <- c_liver_tissue / (c_venous_plasma + cguard)
    Kp_pancreas    <- c_pancreas_tissue / (c_venous_plasma + cguard)
    Kp_spleen      <- c_spleen_tissue / (c_venous_plasma + cguard)
    Kp_gut         <- c_gut_tissue / (c_venous_plasma + cguard)
    Kp_lung        <- c_lung_tissue / (c_venous_plasma + cguard)

    # Eq 40: volume of distribution (L/kg) - total amount in blood and
    # tissues divided by the venous plasma concentration and body weight.
    # The portal vein and the gut lumen are excluded, exactly as in Eq 40.
    a_body <- venous_plasma + venous_bc + arterial_plasma + arterial_bc +
      adipose_iw + adipose_ew + adipose_plasma + adipose_bc +
      bone_iw + bone_ew + bone_plasma + bone_bc +
      brain_iw + brain_ew + brain_plasma + brain_bc +
      heart_iw + heart_ew + heart_plasma + heart_bc +
      kidney_iw + kidney_ew + kidney_plasma + kidney_bc +
      muscle_iw + muscle_ew + muscle_plasma + muscle_bc +
      skin_iw + skin_ew + skin_plasma + skin_bc +
      liver_iw + liver_ew + liver_plasma + liver_bc +
      pancreas_iw + pancreas_ew + pancreas_plasma + pancreas_bc +
      spleen_iw + spleen_ew + spleen_plasma + spleen_bc +
      gut_iw + gut_ew + gut_plasma + gut_bc +
      lung_iw + lung_ew + lung_plasma + lung_bc
    Vdt <- a_body / ((c_venous_plasma + cguard) * bw)

    # ============== 7. Observation ==============
    # Venous plasma concentration. The paper reports no residual-error
    # model: this is a deterministic forward-simulation model.
    Cc <- c_venous_plasma
  })
}
