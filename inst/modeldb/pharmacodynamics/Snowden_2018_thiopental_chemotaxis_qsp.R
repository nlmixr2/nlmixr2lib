Snowden_2018_thiopental_chemotaxis_qsp <- function() {
  description <- paste(
    "QSP. Methodology reference (linked whole-body PBPK and bacterial",
    "chemotaxis signalling model; 25 ODEs). The unreduced linked model of",
    "Snowden, van der Graaf and Tindall 2018, a paper on model reduction,",
    "in which a hypothetical acid given orally is parameterised as",
    "thiopental and acts as the extracellular attractant of the E. coli",
    "CheA/CheY/CheZ/CheB phosphorelay model of Tindall et al. 2009. The",
    "PBPK half is the Jones and Rowland-Yeo perfusion-limited framework",
    "as the authors ran it: 14 states (oral dosing compartment, eleven",
    "tissues, venous and arterial blood), with no testes or rest-of-body",
    "compartment, so 11.5% of the arterial outflow has no return path,",
    "and with hepatic metabolic clearance effectively zero. Both features",
    "are reproduced deliberately because they are what generated the",
    "published Figures 5 and 6. The chemotaxis half is the full",
    "11-species mass-action model; the drug concentration in the gut",
    "tissue drives receptor occupancy through a Hill term, and the output",
    "is the total phosphorylated CheY. RECORD OF THE PUBLISHED SIMULATION,",
    "NOT A USABLE PHARMACOKINETIC MODEL OF thiopental: the compound is",
    "hypothetical, the arterial leak and the missing hepatic clearance",
    "are not physiology, and the PK concentration (mg/L) enters the",
    "ligand Hill term as if it were in umol/L. Deterministic: no",
    "between-subject variability and no residual error are reported.",
    sep = " "
  )
  reference <- paste(
    "Snowden TJ, van der Graaf PH, Tindall MJ. Model reduction in",
    "mathematical pharmacology: integration, reduction and linking of PBPK",
    "and systems biology models. J Pharmacokinet Pharmacodyn.",
    "2018;45(4):537-555. doi:10.1007/s10928-018-9584-y. PMID 29582349.",
    "PMCID PMC6061126. Physiology from main-text Table 1, compound",
    "parameters from Table 2, PBPK equations from Electronic",
    "Supplementary Material 1 (ESM1) Eq. 15 and the printed midazolam",
    "system matrix of ESM1 section 2.1; chemotaxis reactions, equations,",
    "parameters and initial conditions from ESM1 Eqs. 20-23 and ESM1",
    "Table 1, which reproduce Tindall MJ, Porter SL, Wadhams GH, Maini",
    "PK, Armitage JP. Prog Biophys Mol Biol. 2009;100(1-3):40-46.",
    sep = " "
  )
  vignette <- "Snowden_2018_pbpk_chemotaxis_reduction"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  dosing <- c("depot", "venous")

  paper_specific_compartments <- c(
    "chea",
    "cheap",
    "cheap_chey",
    "chea_cheyp",
    "chea_chey",
    "cheyp_chez",
    "chey",
    "cheyp",
    "chez",
    "cheb",
    "chebp"
  )

  compartmentData <- list(
    depot = list(
      analyte = "hypothetical acid (thiopental parameterisation)",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    adipose = list(
      analyte = "hypothetical acid (thiopental parameterisation)",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    ),
    bone = list(
      analyte = "hypothetical acid (thiopental parameterisation)",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    ),
    brain = list(
      analyte = "hypothetical acid (thiopental parameterisation)",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    ),
    gut = list(
      analyte = "hypothetical acid (thiopental parameterisation)",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    ),
    heart = list(
      analyte = "hypothetical acid (thiopental parameterisation)",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    ),
    kidney = list(
      analyte = "hypothetical acid (thiopental parameterisation)",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    ),
    liver = list(
      analyte = "hypothetical acid (thiopental parameterisation)",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    ),
    lung = list(
      analyte = "hypothetical acid (thiopental parameterisation)",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    ),
    muscle = list(
      analyte = "hypothetical acid (thiopental parameterisation)",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    ),
    skin = list(
      analyte = "hypothetical acid (thiopental parameterisation)",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    ),
    spleen = list(
      analyte = "hypothetical acid (thiopental parameterisation)",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    ),
    venous = list(
      analyte = "hypothetical acid (thiopental parameterisation)",
      units = "mg",
      specimen = "whole blood",
      verified = TRUE
    ),
    arterial = list(
      analyte = "hypothetical acid (thiopental parameterisation)",
      units = "mg",
      specimen = "whole blood",
      verified = TRUE
    ),
    chea = list(analyte = "CheA (E. coli cytoplasm)", units = "umol/L", specimen = "not applicable", verified = TRUE),
    cheap = list(
      analyte = "CheA-P (E. coli cytoplasm)",
      units = "umol/L",
      specimen = "not applicable",
      verified = TRUE
    ),
    cheap_chey = list(
      analyte = "CheA-P.CheY complex (E. coli cytoplasm)",
      units = "umol/L",
      specimen = "not applicable",
      verified = TRUE
    ),
    chea_cheyp = list(
      analyte = "CheA.CheY-P complex (E. coli cytoplasm)",
      units = "umol/L",
      specimen = "not applicable",
      verified = TRUE
    ),
    chea_chey = list(
      analyte = "CheA.CheY complex (E. coli cytoplasm)",
      units = "umol/L",
      specimen = "not applicable",
      verified = TRUE
    ),
    cheyp_chez = list(
      analyte = "CheY-P.CheZ complex (E. coli cytoplasm)",
      units = "umol/L",
      specimen = "not applicable",
      verified = TRUE
    ),
    chey = list(analyte = "CheY (E. coli cytoplasm)", units = "umol/L", specimen = "not applicable", verified = TRUE),
    cheyp = list(
      analyte = "CheY-P (E. coli cytoplasm)",
      units = "umol/L",
      specimen = "not applicable",
      verified = TRUE
    ),
    chez = list(analyte = "CheZ (E. coli cytoplasm)", units = "umol/L", specimen = "not applicable", verified = TRUE),
    cheb = list(analyte = "CheB (E. coli cytoplasm)", units = "umol/L", specimen = "not applicable", verified = TRUE),
    chebp = list(analyte = "CheB-P (E. coli cytoplasm)", units = "umol/L", specimen = "not applicable", verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species = "None (methodology paper; a hypothetical compound in a 70 kg reference human physiology linked to an E. coli signalling network; no patients and no fitted estimates).",
    n_subjects = 0L,
    n_studies = 0L,
    disease_state = "N/A (deterministic forward simulation; not a fit of any data).",
    dose_range = "Figures 5 and 6 and main-text Table 3: a single 500 mg oral dose. The main text describes the Figure 6 dose as 150 mg; the published Figure 6 curves are reproduced only by 500 mg (see the vignette).",
    regions = "N/A",
    scope_note = paste(
      "Filed under inst/modeldb/pharmacodynamics/ because the compound is",
      "hypothetical: the paper states that its examples 'do not relate to",
      "specific clinical cases of experimentally tested drugs'. The drug",
      "name in the file stem names the parameter set the authors borrowed",
      "(Pilari and Huisinga 2010, with Rodgers-Rowland partition",
      "coefficients), not a validated model of that drug.",
      sep = " "
    ),
    notes = paste(
      "Printed values that the as-run model does not consume, recorded here",
      "so the source is complete: Table 1 MPPGL = 45 mg/g, FVpl = 0.0424,",
      "FVrb = 0.0347, FVte = 0.01 and FVre = 0.0998 L/kg (the testes and",
      "rest-of-body compartments are absent from the as-run model and their",
      "blood flows and partition coefficients are not printed); Table 2",
      "fup = 0.18 and hepatic blood clearance CLbl = 2.02 mL/min/kg (see",
      "clint below); ESM1 Table 1 totals AT = 7.9, BT = 0.28, YT = 9.7 and",
      "ZT = 3.8 umol/L, which the printed initial conditions satisfy",
      "exactly. Renal clearance and the oral fraction F appear in ESM1",
      "Eq. 15 but have no printed value; the printed midazolam system",
      "matrix implies zero renal clearance, and the Table 2 absorption",
      "rate constants alone reproduce the figure timings.",
      sep = " "
    )
  )

  ini({
    # ---- Physiology: main-text Table 1 (70 kg male) ----
    bw <- fixed(70); label("Body weight (kg)") # Table 1, row 'BW' = 70
    qc <- fixed(390); label("Cardiac output (L/h)") # Table 1, row 'QC' = 3.90 x 10^2

    fvol_adipose  <- fixed(0.213);  label("Fractional volume, adipose (L/kg body weight)") # Table 1, row 'FVad' = 2.13 x 10^-1
    fvol_bone     <- fixed(0.0856); label("Fractional volume, bone (L/kg body weight)") # Table 1, row 'FVbo' = 8.56 x 10^-2
    fvol_brain    <- fixed(0.02);   label("Fractional volume, brain (L/kg body weight)") # Table 1, row 'FVbr' = 2 x 10^-2
    fvol_gut      <- fixed(0.0171); label("Fractional volume, gut (L/kg body weight)") # Table 1, row 'FVgu' = 1.71 x 10^-2
    fvol_heart    <- fixed(0.0047); label("Fractional volume, heart (L/kg body weight)") # Table 1, row 'FVhe' = 4.7 x 10^-3
    fvol_kidney   <- fixed(0.0044); label("Fractional volume, kidney (L/kg body weight)") # Table 1, row 'FVki' = 4.4 x 10^-3
    fvol_liver    <- fixed(0.021);  label("Fractional volume, liver (L/kg body weight)") # Table 1, row 'FVli' = 2.1 x 10^-2
    fvol_lung     <- fixed(0.0076); label("Fractional volume, lung (L/kg body weight)") # Table 1, row 'FVlu' = 7.6 x 10^-3
    fvol_muscle   <- fixed(0.4);    label("Fractional volume, muscle (L/kg body weight)") # Table 1, row 'FVmu' = 4 x 10^-1
    fvol_skin     <- fixed(0.0371); label("Fractional volume, skin (L/kg body weight)") # Table 1, row 'FVsk' = 3.71 x 10^-2
    fvol_spleen   <- fixed(0.0026); label("Fractional volume, spleen (L/kg body weight)") # Table 1, row 'FVsp' = 2.6 x 10^-3
    fvol_venous   <- fixed(0.0514); label("Fractional volume, venous blood (L/kg body weight)") # Table 1, row 'FVve' = 5.14 x 10^-2
    fvol_arterial <- fixed(0.0257); label("Fractional volume, arterial blood (L/kg body weight)") # Table 1, row 'FVar' = 2.57 x 10^-2

    fq_adipose <- fixed(0.05);   label("Fractional blood flow, adipose (fraction of cardiac output)") # Table 1, row 'FQad' = 5 x 10^-2
    fq_bone    <- fixed(0.05);   label("Fractional blood flow, bone (fraction of cardiac output)") # Table 1, row 'FQbo' = 5 x 10^-2
    fq_brain   <- fixed(0.12);   label("Fractional blood flow, brain (fraction of cardiac output)") # Table 1, row 'FQbr' = 1.2 x 10^-1
    fq_gut     <- fixed(0.146);  label("Fractional blood flow, gut (fraction of cardiac output)") # Table 1, row 'FQgu' = 1.46 x 10^-1
    fq_heart   <- fixed(0.04);   label("Fractional blood flow, heart (fraction of cardiac output)") # Table 1, row 'FQhe' = 4 x 10^-2
    fq_kidney  <- fixed(0.19);   label("Fractional blood flow, kidney (fraction of cardiac output)") # Table 1, row 'FQki' = 1.9 x 10^-1
    fq_hepatic <- fixed(0.215);  label("Fractional blood flow, total hepatic venous outflow (fraction of cardiac output)") # Table 1, row 'FQh' = 2.15 x 10^-1
    fq_lung    <- fixed(1);      label("Fractional blood flow, lung (fraction of cardiac output)") # Table 1, row 'FQlu' = 1
    fq_muscle  <- fixed(0.17);   label("Fractional blood flow, muscle (fraction of cardiac output)") # Table 1, row 'FQmu' = 1.7 x 10^-1
    fq_skin    <- fixed(0.05);   label("Fractional blood flow, skin (fraction of cardiac output)") # Table 1, row 'FQsk' = 5 x 10^-2
    fq_spleen  <- fixed(0.0172); label("Fractional blood flow, spleen (fraction of cardiac output)") # Table 1, row 'FQsp' = 1.72 x 10^-2

    # ---- Compound: main-text Table 2, Thiopental column (acid) ----
    bp  <- fixed(0.88);       label("Blood-to-plasma concentration ratio") # Table 2, row 'BP', Thiopental column = 0.88
    lka <- fixed(log(5.64)); label("First-order absorption rate constant (log 1/h)") # Table 2, row 'Ka', Thiopental column = 5.64 1/h
    kp_adipose  <- fixed(12.17); label("Tissue-to-plasma partition coefficient, adipose") # Table 2, row 'Kpad', Thiopental column = 12.17
    kp_bone     <- fixed(1.64); label("Tissue-to-plasma partition coefficient, bone") # Table 2, row 'Kpbo', Thiopental column = 1.64
    kp_brain    <- fixed(1.09); label("Tissue-to-plasma partition coefficient, brain") # Table 2, row 'Kpbr', Thiopental column = 1.09
    kp_gut      <- fixed(2.03); label("Tissue-to-plasma partition coefficient, gut") # Table 2, row 'Kpgu', Thiopental column = 2.03
    kp_heart    <- fixed(1.72); label("Tissue-to-plasma partition coefficient, heart") # Table 2, row 'Kphe', Thiopental column = 1.72
    kp_kidney   <- fixed(4.85); label("Tissue-to-plasma partition coefficient, kidney") # Table 2, row 'Kpki', Thiopental column = 4.85
    kp_liver    <- fixed(3.60); label("Tissue-to-plasma partition coefficient, liver") # Table 2, row 'Kpli', Thiopental column = 3.60
    kp_lung     <- fixed(1.72); label("Tissue-to-plasma partition coefficient, lung") # Table 2, row 'Kplu', Thiopental column = 1.72
    kp_muscle   <- fixed(0.78); label("Tissue-to-plasma partition coefficient, muscle") # Table 2, row 'Kpmu', Thiopental column = 0.78
    kp_skin     <- fixed(1.25); label("Tissue-to-plasma partition coefficient, skin") # Table 2, row 'Kpsk', Thiopental column = 1.25
    kp_spleen   <- fixed(0.94); label("Tissue-to-plasma partition coefficient, spleen") # Table 2, row 'Kpsp', Thiopental column = 0.94
    # Hepatic metabolic clearance acting on the liver tissue concentration
    # (ESM1 Eq. 15g term CLInt). Table 2 prints a hepatic blood clearance
    # CLbl = 2.02 mL/min/kg, but the as-run model carried effectively none:
    # the printed midazolam system matrix (ESM1 section 2.1) has liver
    # diagonal -10.96 against a liver-to-vein entry of 10.95, and Figures 5
    # and 6 are reproduced with zero and not with CLbl converted to L/h.
    clint <- fixed(0); label("As-run hepatic metabolic clearance on liver tissue concentration (L/h)") # ESM1 section 2.1 system matrix, row 7 minus row 12 liver entry = 0.01 (rounding)

    # ---- Chemotaxis signalling: ESM1 Table 1 (from Tindall et al. 2009) ----
    k1   <- fixed(3.75);    label("CheA autophosphorylation rate at zero ligand (1/s)") # ESM1 Table 1, row 'k1' = 3.75 s^-1
    k2   <- fixed(2.50e6);  label("CheA-P + CheY association rate (1/(M*s))") # ESM1 Table 1, row 'k2' = 2.50 x 10^6 (M s)^-1
    km2  <- fixed(15);      label("CheA-P.CheY dissociation rate (1/s)") # ESM1 Table 1, row 'k-2' = 15 s^-1
    k3   <- fixed(1.5e7);   label("CheA-P to CheB phosphotransfer rate (1/(M*s))") # ESM1 Table 1, row 'k3' = 1.5 x 10^7 (M s)^-1
    k5   <- fixed(0.7);     label("CheB-P autodephosphorylation rate (1/s)") # ESM1 Table 1, row 'k5' = 0.7 s^-1
    k6   <- fixed(8.5e-2);  label("CheY-P autodephosphorylation rate (1/s)") # ESM1 Table 1, row 'k6' = 8.5 x 10^-2 s^-1
    k7   <- fixed(650);     label("CheA-P.CheY to CheA.CheY-P phosphotransfer rate (1/s)") # ESM1 Table 1, row 'k7' = 650 s^-1
    km7  <- fixed(50);      label("Reverse phosphotransfer rate CheA.CheY-P to CheA-P.CheY (1/s)") # ESM1 Table 1, row 'k-7' = 50 s^-1
    k8   <- fixed(250);     label("CheA.CheY-P dissociation rate (1/s)") # ESM1 Table 1, row 'k8' = 250 s^-1
    km8  <- fixed(2.08e7);  label("CheA + CheY-P association rate (1/(M*s))") # ESM1 Table 1, row 'k-8' = 2.08 x 10^7 (M s)^-1
    k9   <- fixed(7.50e6);  label("CheA + CheY association rate (1/(M*s))") # ESM1 Table 1, row 'k9' = 7.50 x 10^6 (M s)^-1
    km9  <- fixed(15);      label("CheA.CheY dissociation rate (1/s)") # ESM1 Table 1, row 'k-9' = 15 s^-1
    k10  <- fixed(3.75);    label("CheA.CheY to CheA-P.CheY phosphorylation rate at zero ligand (1/s)") # ESM1 Table 1, row 'k10' = 3.75 s^-1
    k11  <- fixed(5.60e6);  label("CheY-P + CheZ association rate (1/(M*s))") # ESM1 Table 1, row 'k11' = 5.60 x 10^6 (M s)^-1
    km11 <- fixed(0.04);    label("CheY-P.CheZ dissociation rate (1/s)") # ESM1 Table 1, row 'k-11' = 0.04 s^-1
    kcat_chez <- fixed(4.90); label("CheZ-catalysed CheY-P dephosphorylation rate (1/s)") # ESM1 Table 1, row 'k12' = 4.90 s^-1
    lkm   <- fixed(log(2.6)); label("Ligand concentration at half-maximal receptor inhibition (log umol/L)") # ESM1 Table 1, row 'K' = 2.6 uM
    lhill <- fixed(log(1.2)); label("Hill coefficient of ligand-receptor inhibition (log)") # ESM1 Table 1, row 'h' = 1.2

    # Initial conditions: the zero-ligand steady state printed under ESM1 Eq. 21
    chea0       <- fixed(2.579); label("Initial CheA (umol/L)") # ESM1 below Eq. 21, x1(0) = 2.579 uM
    cheap0      <- fixed(3.084); label("Initial CheA-P (umol/L)") # ESM1 below Eq. 21, x2(0) = 3.084 uM
    cheap_chey0 <- fixed(0.079); label("Initial CheA-P.CheY (umol/L)") # ESM1 below Eq. 21, x3(0) = 0.079 uM
    chea_cheyp0 <- fixed(0.728); label("Initial CheA.CheY-P (umol/L)") # ESM1 below Eq. 21, x4(0) = 0.728 uM
    chea_chey0  <- fixed(1.430); label("Initial CheA.CheY (umol/L)") # ESM1 below Eq. 21, x5(0) = 1.430 uM
    cheyp_chez0 <- fixed(2.962); label("Initial CheY-P.CheZ (umol/L)") # ESM1 below Eq. 21, x6(0) = 2.962 uM
    chey0       <- fixed(1.383); label("Initial CheY (umol/L)") # ESM1 below Eq. 21, x7(0) = 1.383 uM
    cheyp0      <- fixed(3.118); label("Initial CheY-P (umol/L)") # ESM1 below Eq. 21, x8(0) = 3.118 uM
    chez0       <- fixed(0.838); label("Initial CheZ (umol/L)") # ESM1 below Eq. 21, x9(0) = 0.838 uM
    cheb0       <- fixed(0.004); label("Initial CheB (umol/L)") # ESM1 below Eq. 21, x10(0) = 0.004 uM
    chebp0      <- fixed(0.276); label("Initial CheB-P (umol/L)") # ESM1 below Eq. 21, x11(0) = 0.276 uM
  })

  model({
    ka <- exp(lka)
    km <- exp(lkm)
    hill <- exp(lhill)

    # ================= PBPK physiology (Table 1) =================
    v_adipose  <- bw * fvol_adipose
    v_bone     <- bw * fvol_bone
    v_brain    <- bw * fvol_brain
    v_gut      <- bw * fvol_gut
    v_heart    <- bw * fvol_heart
    v_kidney   <- bw * fvol_kidney
    v_liver    <- bw * fvol_liver
    v_lung     <- bw * fvol_lung
    v_muscle   <- bw * fvol_muscle
    v_skin     <- bw * fvol_skin
    v_spleen   <- bw * fvol_spleen
    v_venous   <- bw * fvol_venous
    v_arterial <- bw * fvol_arterial

    q_adipose <- qc * fq_adipose
    q_bone    <- qc * fq_bone
    q_brain   <- qc * fq_brain
    q_gut     <- qc * fq_gut
    q_heart   <- qc * fq_heart
    q_kidney  <- qc * fq_kidney
    q_hepatic <- qc * fq_hepatic
    q_lung    <- qc * fq_lung
    q_muscle  <- qc * fq_muscle
    q_skin    <- qc * fq_skin
    q_spleen  <- qc * fq_spleen
    # Hepatic artery: total hepatic outflow less the portal inflow from gut
    # and spleen (ESM1 Eq. 15g QHepArt; matrix entry 11.21 = q_ha / v_arterial).
    q_ha <- q_hepatic - q_gut - q_spleen

    c_adipose  <- adipose / v_adipose
    c_bone     <- bone / v_bone
    c_brain    <- brain / v_brain
    c_gut      <- gut / v_gut
    c_heart    <- heart / v_heart
    c_kidney   <- kidney / v_kidney
    c_liver    <- liver / v_liver
    c_lung     <- lung / v_lung
    c_muscle   <- muscle / v_muscle
    c_skin     <- skin / v_skin
    c_spleen   <- spleen / v_spleen
    c_venous   <- venous / v_venous
    c_arterial <- arterial / v_arterial

    # Venous return from the eight systemic tissues and the liver (ESM1
    # Eq. 15m without its testes and rest-of-body terms, which the as-run
    # model did not carry). The lung still receives the whole cardiac
    # output, so 11.5% of the arterial outflow is lost: this is the as-run
    # structure, visible as a column sum of -24.84 in the printed matrix.
    venous_return <-
      q_adipose * bp * c_adipose / kp_adipose +
      q_bone * bp * c_bone / kp_bone +
      q_brain * bp * c_brain / kp_brain +
      q_heart * bp * c_heart / kp_heart +
      q_kidney * bp * c_kidney / kp_kidney +
      q_hepatic * bp * c_liver / kp_liver +
      q_muscle * bp * c_muscle / kp_muscle +
      q_skin * bp * c_skin / kp_skin

    # ================= PBPK mass balance (mg/h) =================
    # ESM1 Eq. 15 written in amounts; every coefficient matches the printed
    # midazolam system matrix of ESM1 section 2.1, which is in amount form.
    d/dt(depot)    <- -ka * depot                                                                 # ESM1 Eq. 15p (oral dosing compartment)
    d/dt(adipose)  <- q_adipose * c_arterial - q_adipose * bp * c_adipose / kp_adipose           # ESM1 Eq. 15a
    d/dt(bone)     <- q_bone * c_arterial - q_bone * bp * c_bone / kp_bone                       # ESM1 Eq. 15b
    d/dt(brain)    <- q_brain * c_arterial - q_brain * bp * c_brain / kp_brain                   # ESM1 Eq. 15c
    d/dt(gut)      <- q_gut * c_arterial - q_gut * bp * c_gut / kp_gut + ka * depot              # ESM1 Eq. 15d
    d/dt(heart)    <- q_heart * c_arterial - q_heart * bp * c_heart / kp_heart                   # ESM1 Eq. 15e
    d/dt(kidney)   <- q_kidney * c_arterial - q_kidney * bp * c_kidney / kp_kidney               # ESM1 Eq. 15f (renal clearance unprinted; zero in the printed matrix)
    d/dt(liver)    <- q_ha * c_arterial + q_gut * bp * c_gut / kp_gut +
      q_spleen * bp * c_spleen / kp_spleen -
      q_hepatic * bp * c_liver / kp_liver - clint * c_liver                                      # ESM1 Eq. 15g
    d/dt(lung)     <- q_lung * c_venous - q_lung * bp * c_lung / kp_lung                         # ESM1 Eq. 15h
    d/dt(muscle)   <- q_muscle * c_arterial - q_muscle * bp * c_muscle / kp_muscle               # ESM1 Eq. 15i
    d/dt(skin)     <- q_skin * c_arterial - q_skin * bp * c_skin / kp_skin                       # ESM1 Eq. 15j
    d/dt(spleen)   <- q_spleen * c_arterial - q_spleen * bp * c_spleen / kp_spleen               # ESM1 Eq. 15k
    d/dt(venous)   <- venous_return - q_lung * c_venous                                          # ESM1 Eq. 15m (IV input u2 enters here)
    d/dt(arterial) <- q_lung * bp * c_lung / kp_lung - q_lung * c_arterial                       # ESM1 Eq. 15n

    # ================= Linkage (ESM1 section 4) =================
    # The extracellular attractant concentration L is the drug concentration
    # in the gut tissue (ESM1 section 4: 'the predicted concentration of the
    # drug in the gut compartment'). The authors fed the mg/L value into the
    # umol/L Hill term unconverted; that is reproduced here.
    ligand <- c_gut
    u_ligand <- 1 - ligand^hill / (km^hill + ligand^hill)                                        # ESM1 Eq. 23

    # ================= Chemotaxis signalling (ESM1 Eq. 21) =================
    # Rates are printed per second and per molar; model time is hours and
    # species are umol/L, so first-order rates are multiplied by 3600 and
    # second-order rates by 3600 * 1e-6.
    s2h <- 3600
    r1   <- k1 * s2h
    r2   <- k2 * 1e-6 * s2h
    rm2  <- km2 * s2h
    r3   <- k3 * 1e-6 * s2h
    r5   <- k5 * s2h
    r6   <- k6 * s2h
    r7   <- k7 * s2h
    rm7  <- km7 * s2h
    r8   <- k8 * s2h
    rm8  <- km8 * 1e-6 * s2h
    r9   <- k9 * 1e-6 * s2h
    rm9  <- km9 * s2h
    r10  <- k10 * s2h
    r11  <- k11 * 1e-6 * s2h
    rm11 <- km11 * s2h
    r12  <- kcat_chez * s2h

    chea(0)       <- chea0
    cheap(0)      <- cheap0
    cheap_chey(0) <- cheap_chey0
    chea_cheyp(0) <- chea_cheyp0
    chea_chey(0)  <- chea_chey0
    cheyp_chez(0) <- cheyp_chez0
    chey(0)       <- chey0
    cheyp(0)      <- cheyp0
    chez(0)       <- chez0
    cheb(0)       <- cheb0
    chebp(0)      <- chebp0

    # ESM1 Eqs. 21c and 21e print the ligand-modulated CheA.CheY
    # phosphorylation as k1 u; ESM1 Eqs. 20c and 22 name it k10(L) = k1 u,
    # and ESM1 Table 1 gives k10 = k1 = 3.75 s^-1, so k10 is used.
    d/dt(chea)       <- -r1 * u_ligand * chea + r3 * cheap * cheb + r8 * chea_cheyp -
      rm8 * chea * cheyp - r9 * chea * chey + rm9 * chea_chey                                    # ESM1 Eq. 21a
    d/dt(cheap)      <- r1 * u_ligand * chea - r2 * cheap * chey + rm2 * cheap_chey -
      r3 * cheap * cheb                                                                          # ESM1 Eq. 21b
    d/dt(cheap_chey) <- r2 * cheap * chey - rm2 * cheap_chey - r7 * cheap_chey +
      rm7 * chea_cheyp + r10 * u_ligand * chea_chey                                              # ESM1 Eq. 21c
    d/dt(chea_cheyp) <- -r6 * chea_cheyp + r7 * cheap_chey - rm7 * chea_cheyp -
      r8 * chea_cheyp + rm8 * chea * cheyp                                                       # ESM1 Eq. 21d
    d/dt(chea_chey)  <- r6 * chea_cheyp + r9 * chea * chey - rm9 * chea_chey -
      r10 * u_ligand * chea_chey                                                                 # ESM1 Eq. 21e
    d/dt(cheyp_chez) <- r11 * cheyp * chez - rm11 * cheyp_chez - r12 * cheyp_chez                # ESM1 Eq. 21f
    d/dt(chey)       <- -r2 * cheap * chey + rm2 * cheap_chey - r9 * chea * chey +
      rm9 * chea_chey + r12 * cheyp_chez + r6 * cheyp                                            # ESM1 Eq. 21g
    d/dt(cheyp)      <- r8 * chea_cheyp - rm8 * chea * cheyp - r11 * cheyp * chez +
      rm11 * cheyp_chez - r6 * cheyp                                                             # ESM1 Eq. 21h
    d/dt(chez)       <- -r11 * cheyp * chez + rm11 * cheyp_chez + r12 * cheyp_chez               # ESM1 Eq. 21i
    d/dt(cheb)       <- -r3 * cheap * cheb + r5 * chebp                                          # ESM1 Eq. 21j
    d/dt(chebp)      <- r3 * cheap * cheb - r5 * chebp                                           # ESM1 Eq. 21k

    # ================= Outputs =================
    # Cc is the quantity Figure 5 plots as 'plasma concentration': the venous
    # blood concentration multiplied (not divided) by BP. That reading is the
    # only one that reproduces all three Figure 5 panels; see the vignette.
    c_blood <- c_venous
    Cc <- c_venous * bp
    # Total phosphorylated CheY, free and complexed (Figure 6 y-axis).
    cheyp_total <- chea_cheyp + cheyp_chez + cheyp
  })
}
