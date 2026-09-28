Zhou_2021_enrofloxacin_pig_pbpk <- function() {
  description <- paste(
    "Veterinary (pig). PBPK (whole-body, flow-limited, acslXtreme 3.0) for",
    "enrofloxacin given orally to swine as an amorphous solid-dispersion",
    "granule, built to predict the drug concentration in the small-intestinal",
    "contents (the site of action against enteric Campylobacter, Salmonella",
    "and E. coli) and to set the dose (Zhou et al. 2021, Pharmaceutics",
    "13:602). Oral drug empties from the stomach into the small-intestinal",
    "contents by first-order gastric emptying; from there it is absorbed",
    "straight into the liver in competition with first-order faecal loss.",
    "Liver, kidney, muscle, fat and a lumped rest of body are perfused in",
    "parallel from arterial blood and drain to venous blood; the lung sits",
    "in series between venous and arterial blood and receives the whole",
    "cardiac output. Only the unbound fraction of arterial drug perfuses",
    "the tissues. Elimination is first-order hepatic metabolism of the",
    "total liver amount (rate constant proportional to body weight) plus",
    "urinary excretion from the kidney. An intravenous dose goes into",
    "venous blood. Deterministic structure with no between-animal",
    "variability and no residual-error model (the propSd terms are",
    "placeholders). The structure and all but two values are inherited",
    "from the Lin 2016 swine enrofloxacin PBPK code; Zhou 2021 re-calibrated",
    "the gastric-emptying rate constant and the small-intestine volume by",
    "hand. The published intestinal predictions are not reproduced by the",
    "printed absorption rate constant; see the vignette Errata."
  )
  reference <- paste(
    "Zhou K, Huo M, Ma W, Mi K, Xu X, Algharib SA, Xie S, Huang L.",
    "Application of a physiologically based pharmacokinetic model to develop",
    "a veterinary amorphous enrofloxacin solid dispersion. Pharmaceutics.",
    "2021;13(5):602. doi:10.3390/pharmaceutics13050602.",
    "Equations transcribed from the acslXtreme code in the Supplementary",
    "Materials; parameter values from Supplementary Tables S1 and S2 and the",
    "same code. The code and most parameter values are taken by the authors",
    "from Lin Z, Vahl CI, Riviere JE. Human food safety implications of",
    "variation in food animal drug metabolism. Sci Rep. 2016;6:27907.",
    sep = " "
  )
  vignette <- "Zhou_2021_enrofloxacin_pig_pbpk"

  # Doses in mg (the paper prescribes mg/kg; multiply by WT). States hold
  # amounts in mg and volumes are in L (tissue density 1 kg/L), so
  # amount / volume is mg/L = ug/mL = ug/g, the paper's reporting units.
  # The acslX code integrates amounts in umol instead, converting the dose
  # with `MWmol = 2.78` umol/mg and the outputs back with `MWmg = 0.36`
  # mg/umol; the product of those two rounded constants is 1.0008, so every
  # concentration of the published code is 0.08% above the mg-based value
  # here (the exact product 1000/359.4 * 0.3594 is 1).
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")
  # Oral doses go into `stomach` and intravenous doses into `a_venous`;
  # neither is `depot` or `central`.
  dosing <- c("stomach", "a_venous")

  compartmentData <- list(
    # Luminal drug: gastric contents and small-intestinal contents. The
    # small-intestinal contents were sampled in Zhou 2021 (Figure 8), but
    # hold unabsorbed dose, so the closest specimen term is the
    # administration site.
    stomach = list(analyte = "enrofloxacin", units = "mg", specimen = "administration site", verified = TRUE),
    depot = list(analyte = "enrofloxacin", units = "mg", specimen = "administration site", verified = TRUE),
    a_liver = list(analyte = "enrofloxacin", units = "mg", specimen = "tissue", verified = TRUE),
    a_kidney = list(analyte = "enrofloxacin", units = "mg", specimen = "tissue", verified = TRUE),
    a_muscle = list(analyte = "enrofloxacin", units = "mg", specimen = "tissue", verified = TRUE),
    a_fat = list(analyte = "enrofloxacin", units = "mg", specimen = "tissue", verified = TRUE),
    a_remainder = list(analyte = "enrofloxacin", units = "mg", specimen = "tissue", verified = TRUE),
    a_pulmonary = list(analyte = "enrofloxacin", units = "mg", specimen = "tissue", verified = TRUE),
    # The code calls venous and arterial blood "blood/plasma" and compares
    # the venous concentration with the measured plasma concentration.
    a_venous = list(analyte = "enrofloxacin", units = "mg", specimen = "plasma", verified = TRUE),
    a_arterial = list(analyte = "enrofloxacin", units = "mg", specimen = "plasma", verified = TRUE),
    a_urine = list(analyte = "enrofloxacin", units = "mg", specimen = "urine", verified = TRUE),
    a_metabolized = list(analyte = "enrofloxacin", units = "mg", specimen = "not applicable", verified = TRUE),
    a_feces = list(analyte = "enrofloxacin", units = "mg", specimen = "faeces", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      reference_value = 55,
      notes = paste(
        "Every tissue volume and blood flow is a fixed fraction of body",
        "weight and the urinary clearance is per kg, which alone would make",
        "concentrations after a per-kg dose independent of weight. The",
        "hepatic metabolic rate constant, however, is itself proportional",
        "to body weight (code `Km = KmC * BW`, 1/h), so hepatic clearance",
        "grows as WT^2 and heavier pigs have lower plasma exposure for the",
        "same mg/kg dose. The small-intestinal contents are not affected by",
        "WT (their concentration, amount / (0.036 * WT), is weight-invariant).",
        "The code sets `BW = 55` kg ('study-specific; the actual value in",
        "present study'), which reproduces the published intestinal",
        "predictions; Supplementary Table S1 prints a model value of 20 kg",
        "and the pharmacokinetic pigs weighed 25 +/- 5 kg (Zhou 2021",
        "Section 2.2). reference_value is the code default."
      ),
      source_name = "BW"
    )
  )

  population <- list(
    species = "pig (three-way crossbred swine)",
    n_subjects = 21L,
    n_studies = 2L,
    age_range = NA_character_,
    weight_range = "25 +/- 5 kg (Zhou 2021 Section 2.2); the code sets BW = 55 kg for the intestinal-content simulations",
    sex_female_pct = NA_real_,
    disease_state = "healthy",
    dose_range = "Single intragastric 2.5 mg/kg of the enrofloxacin solid-dispersion granule (plasma validation, n = 6); 5 mg/kg twice daily for 5 days (intestinal-content validation, n = 15); 5 and 10 mg/kg once daily and 5 mg/kg twice daily simulated for dose selection",
    regions = "China (Huazhong Agricultural University, Wuhan)",
    notes = paste(
      "Eighteen clinically healthy three-way crossbred pigs were bought for",
      "the study (Section 2.2). Six received a single intragastric dose of",
      "2.5 mg/kg and were bled at 0.25-48 h (Section 2.6, Table 2); the",
      "measured plasma concentrations validated the model (Figure 7). The",
      "intestinal-content samples (1, 109, 112, 120 and 132 h, three pigs",
      "each) came from fifteen pigs of a later tissue-residue study dosed",
      "at 5 mg/kg twice daily for five days (Figure 8). No parameter was",
      "estimated statistically: the structure and values are those of the",
      "Lin 2016 swine enrofloxacin PBPK model (Zhou 2021 ref 38), with the",
      "gastric-emptying rate constant and the small-intestine volume",
      "adjusted by hand to the observed data (Section 2.5). There are no",
      "random effects; the model describes a typical pig."
    )
  )

  ini({
    # ---------------------------------------------------------------
    # Oral absorption (Supplementary Table S2 'Absorption rate constant'
    # and the acslX code). Kst was re-calibrated by hand in Zhou 2021
    # (Section 2.5: 'visually reasonable values for ... gastric emptying
    # rate constant were obtained by an iterative manual adjustment'), so
    # it is left unfixed; Ka and Kfeces are the Lin 2016 values.
    # ---------------------------------------------------------------
    lkst <- log(2)
    label("Gastric emptying rate constant Kst, stomach to small-intestinal contents (1/h)") # Table S2 (Kst model value 2.0; Lin 2016 value 1.0); code `Kst = 2`
    lka <- fixed(log(0.55))
    label("Intestinal absorption rate constant Ka, small-intestinal contents to liver (1/h)") # Table S2 (Ka 0.55); code `Ka = 0.55`. The published intestinal predictions need Ka near 0.25; see vignette Errata
    lkfec <- fixed(log(0.01))
    label("Faecal elimination rate constant Kfeces for unabsorbed drug (1/h)") # Table S2 (Kfeces 0.01); code `Kfeces = 0.01`

    # ---------------------------------------------------------------
    # Tissue:plasma partition coefficients (Table S2 and the code;
    # Buur 2005 values via Lin 2016). Prest is in the code only.
    # ---------------------------------------------------------------
    lkp_liver <- fixed(log(4.3))
    label("Liver:plasma partition coefficient PL (unitless)") # Table S2 (PL 4.3); code `PL = 4.3`
    lkp_kidney <- fixed(log(5.5))
    label("Kidney:plasma partition coefficient PK (unitless)") # Table S2 (PK 5.5); code `PK = 5.5`
    lkp_muscle <- fixed(log(3))
    label("Muscle:plasma partition coefficient PM (unitless)") # Table S2 (PM 3.0); code `PM = 3`
    lkp_fat <- fixed(log(0.53))
    label("Fat:plasma partition coefficient PF (unitless)") # Table S2 (PF 0.53); code `PF = 0.53`
    lkp_lung <- fixed(log(4.3))
    label("Lung:plasma partition coefficient PLu (unitless)") # Table S2 (Plu 4.3); code `PLu = 4.3`
    lkp_remainder <- fixed(log(8))
    label("Rest-of-body:plasma partition coefficient Prest (unitless)") # code `Prest = 8` (not tabulated in Table S2)

    # ---------------------------------------------------------------
    # Protein binding and elimination (Table S2 and the code).
    # ---------------------------------------------------------------
    fu <- fixed(0.54)
    label("Unbound fraction in arterial blood, 1 - PB (unitless)") # Table S2 (PB 0.46 bound); code `PB = 0.46`, `CAfree = CA*(1-PB)`
    lkmet <- fixed(log(0.045))
    label("Hepatic metabolic rate constant per kg body weight KmC (1/h per kg)") # Table S2 (KmC 0.045 /(h*kg)); code `KmC = 0.045`, `Km = KmC*BW`
    lcl_renal <- fixed(log(0.12))
    label("Urinary clearance per kg body weight KurineC (L/h/kg), applied to kidney venous concentration") # Table S2 (KurineC 0.12 L/h/kg); code `KurineC = 0.12`, `Rurine = Kurine*CVK`

    # ---------------------------------------------------------------
    # Zhou 2021 calibrated the model by eye and checked it by linear
    # regression of predicted on observed (Figures 7B, 8B); it reports no
    # residual-error model. nlmixr2 model definitions require a
    # residual-error term, so the propSd terms below are fixed
    # placeholders for syntactic completeness only and must NOT be read as
    # estimates. Same convention as Mi_2023_cefquinome_pbpk.
    # ---------------------------------------------------------------
    propSd <- fixed(0.10)
    label("Proportional residual error placeholder, plasma (fraction)") # not reported in Zhou 2021; placeholder only
    propSd_Cintestine <- fixed(0.10)
    label("Proportional residual error placeholder, small-intestinal contents (fraction)") # not reported in Zhou 2021; placeholder only
  })

  model({
    # =================================================================
    # Swine physiology (Supplementary Table S1 and the acslX INITIAL
    # block; Buur 2005 / Upton 2008 via Lin 2016, small intestine from
    # Lautz 2020 re-calibrated by hand). These are physiology rather than
    # drug parameters, so they are traceable literals here, as in
    # Mi_2023_cefquinome_pbpk.
    # =================================================================
    qcc <- 5 # cardiac output, L/h/kg; Table S1 (QCC 5.0); code `QCC = 5`
    qc_liver <- 0.2725 # fraction of cardiac output; Table S1 (QLC 0.2725); code `QLC = 0.2725`
    qc_kidney <- 0.12 # Table S1 (QKC 0.12); code `QKC = 0.12`
    qc_muscle <- 0.251 # Table S1 (QMC 0.251); code `QMC = 0.251`
    qc_fat <- 0.1275 # Table S1 (QFC 0.1275); code `QFC = 0.1275`

    vc_liver <- 0.0247 # fraction of body weight; Table S1 (VLC 0.0247); code `VLC = 0.0247`
    vc_kidney <- 0.004 # Table S1 (VKC 0.004); code `VKC = 0.004`
    vc_muscle <- 0.4 # Table S1 (VMC 0.4); code `VMC = 0.40`
    vc_fat <- 0.32 # Table S1 (VFC 0.32); code `VFC = 0.32`
    vc_lung <- 0.01 # Table S1 (VLuC 0.01); code `VLuC = 0.01`
    vc_blood <- 0.06 # code `VBloodC = 0.06`; Table S1 splits it as arterial 0.0156 + venous 0.0444
    # Small-intestine contents volume, re-calibrated by hand (Section 2.5);
    # Table S1 (VSiC model value 0.036, published 0.0126); code `VSiC = 0.036`
    vc_small_intestine <- 0.036

    # =================================================================
    # Drug-specific parameters
    # =================================================================
    kst <- exp(lkst)
    ka <- exp(lka)
    kfec <- exp(lkfec)
    kp_liver <- exp(lkp_liver)
    kp_kidney <- exp(lkp_kidney)
    kp_muscle <- exp(lkp_muscle)
    kp_fat <- exp(lkp_fat)
    kp_lung <- exp(lkp_lung)
    kp_remainder <- exp(lkp_remainder)
    kmet <- exp(lkmet) * WT # code `Km = KmC*BW` (1/h)
    cl_renal <- exp(lcl_renal) * WT # code `Kurine = KurineC*BW` (L/h)

    # =================================================================
    # Flows (L/h) and volumes (L). The rest-of-body flow and volume are
    # the complements, code `Qrest = QC-QL-QK-QM-QF` (0.229 of cardiac
    # output, Table S1 QrestC) and `Vrest = BW-VL-VK-VM-VF-VLu-VBlood`
    # (0.1813 of body weight, Table S1 VrestC). The lung carries the whole
    # cardiac output. Blood is split 74% venous, 26% arterial (code
    # `Vven = VBlood*0.74`, `Vart = VBlood*0.26`). The intestinal contents
    # are not part of body weight in the complement.
    # =================================================================
    qc <- qcc * WT
    q_liver <- qc_liver * qc
    q_kidney <- qc_kidney * qc
    q_muscle <- qc_muscle * qc
    q_fat <- qc_fat * qc
    q_remainder <- qc - q_liver - q_kidney - q_muscle - q_fat

    v_liver <- vc_liver * WT
    v_kidney <- vc_kidney * WT
    v_muscle <- vc_muscle * WT
    v_fat <- vc_fat * WT
    v_lung <- vc_lung * WT
    v_blood <- vc_blood * WT
    v_venous <- v_blood * 0.74
    v_arterial <- v_blood * 0.26
    v_remainder <- WT - v_liver - v_kidney - v_muscle - v_fat - v_lung - v_blood
    v_small_intestine <- vc_small_intestine * WT

    # =================================================================
    # Concentrations (mg/L). Tissue outflow is at the venous equilibrium
    # concentration C / P. Arterial drug perfuses the tissues UNBOUND
    # only (code `CAfree = CA*(1-PB)`), while every tissue drains its
    # total venous-equilibrium concentration -- exactly as coded.
    # =================================================================
    c_venous <- a_venous / v_venous
    c_arterial <- a_arterial / v_arterial
    c_arterial_free <- c_arterial * fu
    c_lung <- a_pulmonary / v_lung
    c_liver <- a_liver / v_liver
    c_kidney <- a_kidney / v_kidney
    c_muscle <- a_muscle / v_muscle
    c_fat <- a_fat / v_fat
    c_remainder <- a_remainder / v_remainder
    cv_lung <- c_lung / kp_lung
    cv_liver <- c_liver / kp_liver
    cv_kidney <- c_kidney / kp_kidney
    cv_muscle <- c_muscle / kp_muscle
    cv_fat <- c_fat / kp_fat
    cv_remainder <- c_remainder / kp_remainder

    # =================================================================
    # Fluxes (mg/h)
    # =================================================================
    r_absorb <- ka * depot # code `RAO = Ka*AI`
    r_feces <- kfec * depot # code `Rfeces = Kfeces*AI`
    r_met <- kmet * a_liver # code `Rmet = Km*CL*VL`
    r_renal <- cl_renal * cv_kidney # code `Rurine = Kurine*CVK`

    # =================================================================
    # ODEs (acslX DERIVATIVE block)
    # =================================================================
    d/dt(stomach) <- -kst * stomach # code `RAST = RDOSEoral - Kst*AST`
    d/dt(depot) <- kst * stomach - r_absorb - r_feces # code `RAI = Kst*AST - Ka*AI - Kfeces*AI`
    # Absorbed drug enters the liver directly (no portal or gut tissue).
    d/dt(a_liver) <- q_liver * (c_arterial_free - cv_liver) + r_absorb - r_met # code `RL`
    d/dt(a_kidney) <- q_kidney * (c_arterial_free - cv_kidney) - r_renal # code `RK`
    d/dt(a_muscle) <- q_muscle * (c_arterial_free - cv_muscle) # code `RM`
    d/dt(a_fat) <- q_fat * (c_arterial_free - cv_fat) # code `RF`
    d/dt(a_remainder) <- q_remainder * (c_arterial_free - cv_remainder) # code `Rrest`
    d/dt(a_pulmonary) <- qc * (c_venous - cv_lung) # code `RALu = QC*(CV-CVLu)`
    d/dt(a_venous) <- q_liver * cv_liver + q_kidney * cv_kidney + q_muscle * cv_muscle +
      q_fat * cv_fat + q_remainder * cv_remainder - qc * c_venous # code `RV` (IV dose enters here)
    d/dt(a_arterial) <- qc * cv_lung - qc * c_arterial_free # code `RA = QC*CVLu - QC*CAfree`
    d/dt(a_urine) <- r_renal # code `Aurine = Integ(Rurine)`
    d/dt(a_metabolized) <- r_met # code `Amet = Integ(Rmet)`
    d/dt(a_feces) <- r_feces # code `Afeces = Integ(Rfeces)`

    # =================================================================
    # Outputs (ug/mL = mg/L; ug/g for tissues)
    # =================================================================
    Cc <- c_venous # plasma; code `CVmg = CV*MWmg`
    Cliver <- c_liver # code `CLmg`
    Ckidney <- c_kidney # code `Ckmg`
    Cmuscle <- c_muscle # code `CMmg`
    Cfat <- c_fat # code `CFmg`
    Clung <- c_lung
    # Concentration in the small-intestinal contents, code `CAI = AI/VSi`
    # converted to mg/L. The code's own mg output line `CAImg = AI*MWmg`
    # omits the division by VSi, so the 'intestinal concentration' Zhou
    # 2021 reports (Section 3.7, Figure 8) is numerically the AMOUNT in
    # the small-intestinal contents in mg, i.e. the `depot` state. See the
    # vignette Errata.
    Cintestine <- depot / v_small_intestine

    Cc ~ prop(propSd)
    Cintestine ~ prop(propSd_Cintestine)
  })
}
