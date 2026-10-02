Humphries_2021_kanamycin <- function() {
  description <- paste(
    "PBPK (permeability-limited lung, 25 ODEs; Simcyp V16 R1). Kanamycin",
    "plasma, epithelial lining fluid and lung tissue concentrations after an",
    "intramuscular dose, from the anti-tuberculosis PBPK compound files of",
    "Humphries et al. placed in the multicompartment permeability-limited lung",
    "model of Gaohua et al. 2015 that the paper reuses.  The lung is divided",
    "into seven segments - upper and lower airways plus the five lobes - each",
    "carrying pulmonary capillary blood, tissue mass and fluid (mucus and",
    "epithelial lining fluid, in immediate equilibrium with alveolar air),",
    "plus the pulmonary blood reservoir and arterial blood.  Passive",
    "permeation of unbound unionised drug crosses the apical (fluid-mass) and",
    "basal (mass-blood) membranes, with Henderson-Hasselbalch ionisation at",
    "the fluid (pH 6.6), mass (pH 6.69) and blood (pH 7.4) pH values;",
    "transporters and lung metabolism are exposed and fixed at zero.",
    "DEVIATION FROM THE PUBLISHED ARCHITECTURE: the paper embeds the lung in a",
    "Simcyp full PBPK whose perfusion-limited tissues need per-tissue Kp values",
    "that are not published, so the systemic side is reduced to one",
    "well-stirred compartment at the compound file's own Vss, intramuscular",
    "ka and fa, and renal clearance.  Deterministic typical-value simulation",
    "model: the paper reports no IIV and no residual-error model.  See the",
    "validation vignette for the reduction's accuracy against the paper's own",
    "predicted plasma exposure and lung profile."
  )
  reference <- paste(
    "Humphries H, Almond L, Berg A, Gardner I, Hatley O, Pan X, Small B,",
    "Zhang M, Jamei M, Romero K. Development of physiologically-based",
    "pharmacokinetic models for standard of care and newer tuberculosis drugs.",
    "CPT Pharmacometrics Syst Pharmacol. 2021;10(11):1382-1395.",
    "doi:10.1002/psp4.12707.",
    "Lung model structure and lung physiology: Gaohua L, Wedagedera J, Small",
    "BG, Almond L, Romero K, Hermann D, Hanna D, Jamei M, Gardner I.",
    "Development of a Multicompartment Permeability-Limited Lung PBPK Model and",
    "Its Application in Predicting Pulmonary Pharmacokinetics of",
    "Antituberculosis Drugs. CPT Pharmacometrics Syst Pharmacol.",
    "2015;4(10):605-613. doi:10.1002/psp4.12034 (reference 7 of Humphries",
    "2021), with reference physiology from Jamei M et al. Clin Pharmacokinet.",
    "2014;53:73-87, Electronic Supplementary Material 1.",
    sep = " "
  )
  vignette <- "Humphries_2021_tb_lung_pbpk"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # The lung segments and the pulmonary blood reservoir are anatomical states
  # of the Gaohua 2015 lung model and do not generalise to other extractions;
  # `arterial` is the systemic arterial blood pool of its Appendix S1 eq 7.
  paper_specific_compartments <- c("arterial")
  paper_specific_compartment_pattern <- "^lung_"

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Scales the systemic distribution volume only (Vss is reported in",
        "L/kg in Supplementary Table S2).  Renal clearance and the lung",
        "physiology are absolute values for the reference adult and are not",
        "weight scaled; the paper reports no allometric relationship.  Use 70",
        "kg to reproduce the simulations in this package's vignette.",
        sep = " "
      ),
      source_name = "WT"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "kanamycin", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "kanamycin", units = "mg", specimen = "plasma", verified = TRUE),
    arterial = list(analyte = "kanamycin", units = "mg", specimen = "whole blood", verified = TRUE),
    lung_pbr = list(analyte = "kanamycin", units = "mg", specimen = "whole blood", verified = TRUE),
    lung_rt_fluid = list(analyte = "kanamycin", units = "mg", specimen = "epithelial lining fluid", verified = TRUE),
    lung_rt_mass = list(analyte = "kanamycin", units = "mg", specimen = "tissue", verified = TRUE),
    lung_rt_blood = list(analyte = "kanamycin", units = "mg", specimen = "whole blood", verified = TRUE),
    lung_rm_fluid = list(analyte = "kanamycin", units = "mg", specimen = "epithelial lining fluid", verified = TRUE),
    lung_rm_mass = list(analyte = "kanamycin", units = "mg", specimen = "tissue", verified = TRUE),
    lung_rm_blood = list(analyte = "kanamycin", units = "mg", specimen = "whole blood", verified = TRUE),
    lung_rl_fluid = list(analyte = "kanamycin", units = "mg", specimen = "epithelial lining fluid", verified = TRUE),
    lung_rl_mass = list(analyte = "kanamycin", units = "mg", specimen = "tissue", verified = TRUE),
    lung_rl_blood = list(analyte = "kanamycin", units = "mg", specimen = "whole blood", verified = TRUE),
    lung_lt_fluid = list(analyte = "kanamycin", units = "mg", specimen = "epithelial lining fluid", verified = TRUE),
    lung_lt_mass = list(analyte = "kanamycin", units = "mg", specimen = "tissue", verified = TRUE),
    lung_lt_blood = list(analyte = "kanamycin", units = "mg", specimen = "whole blood", verified = TRUE),
    lung_ll_fluid = list(analyte = "kanamycin", units = "mg", specimen = "epithelial lining fluid", verified = TRUE),
    lung_ll_mass = list(analyte = "kanamycin", units = "mg", specimen = "tissue", verified = TRUE),
    lung_ll_blood = list(analyte = "kanamycin", units = "mg", specimen = "whole blood", verified = TRUE),
    lung_la_fluid = list(analyte = "kanamycin", units = "mg", specimen = "epithelial lining fluid", verified = TRUE),
    lung_la_mass = list(analyte = "kanamycin", units = "mg", specimen = "tissue", verified = TRUE),
    lung_la_blood = list(analyte = "kanamycin", units = "mg", specimen = "whole blood", verified = TRUE),
    lung_ua_fluid = list(analyte = "kanamycin", units = "mg", specimen = "epithelial lining fluid", verified = TRUE),
    lung_ua_mass = list(analyte = "kanamycin", units = "mg", specimen = "tissue", verified = TRUE),
    lung_ua_blood = list(analyte = "kanamycin", units = "mg", specimen = "whole blood", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 390L,
    n_studies = 2L,
    age_range = "21-59 years",
    sex_female_pct = 12.8,
    disease_state = "virtual healthy adults (Simcyp library populations) matched to healthy-volunteer and tuberculosis-patient studies",
    dose_range = "500 mg single intramuscular dose (plasma verification); 1000 mg single dose (lung)",
    regions = "Simcyp Sim-Healthy Volunteer and Sim-North European Caucasian virtual populations",
    notes = paste(
      "Supplementary Table S3: the plasma verification simulation matched",
      "Cabana and Taggart 1973 (24 healthy men, 21-48 years, 500 mg",
      "intramuscular single dose) with 10 trials in the Simcyp Sim-HV",
      "population.  The lung simulation (Figure 5) matched the lesion study",
      "of Prideaux et al. 2015 / Strydom et al. 2019 (15 tuberculosis",
      "patients, 23-59 years, 33% women, 1000 mg single dose) with 10 trials",
      "of 15 Sim-NEC subjects.  These are virtual subjects generated by the",
      "Simcyp population library, not enrolled participants; 390 = 240 + 150",
      "virtual subjects, of whom 50 were women.",
      sep = " "
    )
  )

  ini({
    # All lung-physiology values below are those of the Gaohua 2015 lung
    # model, which Humphries 2021 reuses ("constructed using a previously
    # published multicompartment permeability-limited lung model [ref 7]
    # (Figure S1)"); the source comments therefore point to Gaohua 2015.
    # ================= System: reference physiology =================
    # Cardiac output and the total lung / blood volumes come from the upstream
    # Simcyp full-body PBPK the Appendix defers to (Jamei 2014 ESM1); every
    # regional split below is printed in the main text under "Model
    # parameters".
    qc         <- fixed(356);   label("Cardiac output (L/h)")                              # Gaohua 2015 via Jamei 2014 ESM1 Table S3 ("On average the cardiac output (CO) is 356 (L/h)")
    vLungTotal <- fixed(0.53);  label("Total lung volume (L)")                             # Gaohua 2015 via Jamei 2014 ESM1 Table S1 (Lungs 0.53 L)
    vArterial  <- fixed(1.955); label("Arterial blood volume (L)")                         # Gaohua 2015 via Jamei 2014 ESM1 Table S1 (blood 5.75 L; "venous blood volume is assumed to be 66% of the total blood volume") -> 0.34 * 5.75
    vPbrTotal  <- fixed(0.089); label("Total pulmonary blood reservoir volume (L)")        # Gaohua 2015 Methods "Tissue volumes" (mean 89 mL, CV 20%, literature range 63-114 mL)
    vElfTotal  <- fixed(0.025); label("Total epithelial lining fluid volume (L)")          # Gaohua 2015 Methods "Tissue volumes" (mean 25 mL, CV 20%, literature range 10-40 mL)
    vAlvTotal  <- fixed(5.6);   label("Total alveolar air volume (L)")                     # Gaohua 2015 Methods "Distribution of alveolar volume" (5.6 L, CV 10%)
    vAirwayAir <- fixed(0.05);  label("Air volume of one airway segment (L)")              # Gaohua 2015 Methods "Upper/lower airways" (UA and LA both 0.05 L)

    # PBR, lung mass and ELF are all split on the same proportions: 16.67% to
    # each of the five lobes, the remaining 16.67% shared equally by UA and LA.
    fSegLobe   <- fixed(0.1666667);  label("Fraction of PBR / mass / ELF in each lung lobe")     # Gaohua 2015 Methods "Tissue volumes" (each of the five lobes receives 16.67% of the total PBR)
    fSegAirway <- fixed(0.08333333); label("Fraction of PBR / mass / ELF in each airway segment")# Gaohua 2015 Methods "Tissue volumes" (remaining 16.67% split equally between the UA and LA segments)

    # Alveolar air volume by lobe (Gaohua 2015 Methods "Distribution of alveolar volume").
    fAlvRt <- fixed(0.1578947); label("Fraction of alveolar air volume, right top lobe")    # Gaohua 2015 Methods (3/19 RT)
    fAlvRm <- fixed(0.1052632); label("Fraction of alveolar air volume, right middle lobe") # Gaohua 2015 Methods (2/19 RM)
    fAlvRl <- fixed(0.2631579); label("Fraction of alveolar air volume, right lower lobe")  # Gaohua 2015 Methods (5/19 RL)
    fAlvLt <- fixed(0.2105263); label("Fraction of alveolar air volume, left top lobe")     # Gaohua 2015 Methods (4/19 LT)
    fAlvLl <- fixed(0.2631579); label("Fraction of alveolar air volume, left lower lobe")   # Gaohua 2015 Methods (5/19 LL)

    # Pulmonary blood flow as a fraction of cardiac output (Gaohua 2015 Methods "Blood
    # flow rate"); the upper airway is perfused from arterial blood instead.
    fqRt <- fixed(0.086); label("Pulmonary blood flow fraction, right top lobe")    # Gaohua 2015 Methods "Blood flow rate" (8.6% RT)
    fqRm <- fixed(0.118); label("Pulmonary blood flow fraction, right middle lobe") # Gaohua 2015 Methods "Blood flow rate" (11.8% RM)
    fqRl <- fixed(0.310); label("Pulmonary blood flow fraction, right lower lobe")  # Gaohua 2015 Methods "Blood flow rate" (31.0% RL)
    fqLt <- fixed(0.086); label("Pulmonary blood flow fraction, left top lobe")     # Gaohua 2015 Methods "Blood flow rate" (8.6% LT)
    fqLl <- fixed(0.349); label("Pulmonary blood flow fraction, left lower lobe")   # Gaohua 2015 Methods "Blood flow rate" (34.9% LL)
    fqLa <- fixed(0.050); label("Pulmonary blood flow fraction, lower airway")      # Gaohua 2015 Methods "Blood flow rate" (5.0% LA)
    fqUa <- fixed(0.025); label("Arterial blood flow fraction, upper airway")       # Gaohua 2015 Methods "Blood flow rate" (2.5% of arterial blood circulates in the UA)

    # Alveolar ventilation.  The total is set from the ventilation/perfusion
    # ratio and then partitioned by lobe (Gaohua 2015 Methods "Ventilation/perfusion
    # distribution").
    vqRatio <- fixed(1.0);   label("Ventilation/perfusion ratio (unitless)")   # Gaohua 2015 Methods "Ventilation/perfusion distribution" (geometric mean = 1.0)
    fvRt    <- fixed(0.149); label("Ventilation fraction, right top lobe")     # Gaohua 2015 Methods (14.9% RT)
    fvRm    <- fixed(0.124); label("Ventilation fraction, right middle lobe")  # Gaohua 2015 Methods (12.4% RM)
    fvRl    <- fixed(0.259); label("Ventilation fraction, right lower lobe")   # Gaohua 2015 Methods (25.9% RL)
    fvLt    <- fixed(0.148); label("Ventilation fraction, left top lobe")      # Gaohua 2015 Methods (14.8% LT)
    fvLl    <- fixed(0.320); label("Ventilation fraction, left lower lobe")    # Gaohua 2015 Methods (32.0% LL)

    # Absorption (permeation) surface area: the deep lung is apportioned by
    # relative alveolar volume, the airways are split 50:50 (Methods
    # "Absorption area").
    saDeep   <- fixed(140); label("Deep-lung permeation surface area (m^2)")   # Gaohua 2015 Methods "Absorption area" (~140 m2, CV 30%, assigned to each lobe by relative alveolar volume)
    saAirway <- fixed(1.5); label("Upper + lower airway surface area (m^2)")   # Gaohua 2015 Methods "Absorption area" (1.5 m2, CV 50%, for the sum of the UA and LA, split 50:50)

    # Compartment pH used in the Henderson-Hasselbalch ionisation terms.
    phFluid <- fixed(6.60); label("pH of the lung fluid (ELF) compartment")    # Gaohua 2015 Methods "pH" (pH about 6.6 in healthy individuals)
    phMass  <- fixed(6.69); label("pH of the lung tissue mass compartment")    # Gaohua 2015 Methods "pH" (pH in the lung mass is around 6.69 +/- 0.07)
    phBlood <- fixed(7.40); label("pH of blood / plasma")                      # Gaohua 2015 Methods "pH" (arterial blood pH between 7.38 and 7.43); the ionisation corrections in Results are quoted at pH 7.4

    # Lung metabolism and the four transporter clearances are the paper's base
    # case: passive diffusion only.  They are exposed so a user can reproduce
    # the published sensitivity analyses.
    clMet      <- fixed(0); label("Metabolic clearance in the lung mass compartment (L/h)")            # Gaohua 2015 Methods "Parameterization of the multicompartment lung model" (for all compounds, lung metabolism was assumed to be negligible, CLmet = 0)
    clUptakeFm <- fixed(0); label("Apical (fluid->mass) uptake transporter clearance (L/h)")           # Gaohua 2015 Methods (initial simulations considered distribution by passive diffusion only)
    clEffluxFm <- fixed(0); label("Apical (mass->fluid) efflux transporter clearance (L/h)")           # Gaohua 2015 Methods; Figure 5 varies this over 0, 0.06, 0.6, 6 and 60 L/h
    clUptakeMb <- fixed(0); label("Basal (blood->mass) uptake transporter clearance (L/h)")            # Gaohua 2015 Methods; Supplementary Figure 4 varies this over 0, 0.6 and 60 L/h
    clEffluxMb <- fixed(0); label("Basal (mass->blood) efflux transporter clearance (L/h)")            # Gaohua 2015 Methods (initial simulations considered distribution by passive diffusion only)

    # Gaohua 2015 Appendix S1 requires TWO permeability-surface products per
    # segment -- apical (fluid-mass, CL_PD,FM) and basal (mass-blood,
    # CL_PD,MB) -- but Humphries 2021 Table S2, like Gaohua 2015 Table S2,
    # publishes ONE "lung effective permeability" per compound.  The shipped
    # default of 1 makes the two equal, the reading most consistent with a
    # single published permeability and a single in vitro-in vivo
    # extrapolation.  Exposed because it is an unpublished structural choice.
    ratioPdBasal <- fixed(1); label("Basal:apical permeability-surface ratio (unitless)")  # NOT PUBLISHED: Gaohua 2015 Appendix S1 needs CL_PD,MB and CL_PD,FM separately; Table S2 gives one permeability

    # Air:fluid partition coefficient, K_AF = H / (R * T) (Gaohua 2015
    # Appendix S1 eq 1).  Humphries 2021 Table S2 prints Henry's constant for
    # kanamycin; body temperature is not printed and is taken as 37 degC.
    henry <- fixed(2.95e-33); label("Henry's constant (Pa*m^3/mol)")        # Table S2 "Henry's Constant" (2.95 E-33; footnote d, predicted using EPI Suite)
    tempK <- fixed(310.15);   label("Body temperature (K)")                  # NOT PRINTED: 37 degC assumed; enters only K_AF = H / (R * T), which is ~1e-36 here

    # ================= Compound file: kanamycin =================
    # Humphries 2021 Supplementary Tables S1 and S2 unless noted.  The
    # platform's Kp scalar (0.2, Table S2) acts only through the whole-body
    # Kp values; its effect on the plasma level is carried by the optimised
    # Vss below.
    bp       <- fixed(0.644); label("Blood:plasma concentration ratio")                # Table S1 "BP" (0.644; footnote f, predicted using Simcyp V16)
    fup      <- fixed(0.99);  label("Unbound fraction in plasma")                      # Table S1 "f u,p" (0.99; footnote c, experimental)
    pka1     <- fixed(9.5);   label("pKa (monoprotic base)")                           # Table S1 "pKa" (9.5; footnote d); compound type "Monoprotic Base"
    lka      <- fixed(log(2)); label("First-order intramuscular absorption rate constant (1/h)")  # Table S1 "Absorption Model": "im (venous blood, fa 1, ka 2 h)"
    lfdepot  <- fixed(log(1)); label("Fraction of the intramuscular dose absorbed (unitless)")    # Table S1 "Absorption Model": "im (venous blood, fa 1, ka 2 h)"
    lvc      <- fixed(log(0.236)); label("Systemic distribution volume (L/kg)")        # Table S2 "V SS (L/kg)" (0.236; footnote b, optimized against clinical data)
    lcl_renal <- fixed(log(4.74)); label("Renal clearance (L/h)")                      # Table S2 "CL R (L/h)" (4.74); Table S2 "Elimination" is 0 (no metabolic clearance)

    peffLung <- fixed(0.00617e-4); label("Lung effective permeability (cm/s)")         # Table S2 "Lung effective permeability (10^-4 cm/s)" (0.00617; footnote e, LogD/HBD QSAR corrected for fraction unionised at pH 7.4)
    fuMass   <- fixed(0.999); label("Unbound fraction in the lung tissue mass")        # Table S2 "fu mass" (0.999; footnote f, predicted using Simcyp V16)
    fuFluid  <- fixed(1);     label("Unbound fraction in the lung fluid (ELF)")        # NOT PRINTED in Humphries 2021; Gaohua 2015 Table S2 "fu fluid" = 1 for every compound

    # No residual-error model is reported: the paper is a deterministic
    # forward simulation.  An additive term fixed at zero keeps the model
    # solvable by nlmixr2 without inventing a variance.
    addSd <- fixed(0); label("Additive residual error (mg/L)")                         # not reported; the paper presents predictions, not a fit
  })

  model({
    # =============== 1. Reference-individual scaling ===============
    # The systemic reduction (see description): one well-stirred systemic
    # compartment carrying the whole body except the lung, at the compound
    # file's own Vss, with the lung layer solved as published by Gaohua 2015.
    vc <- exp(lvc) * WT
    ka <- exp(lka)

    # =============== 2. Ionisation, air partition and clearance ===============
    # Kanamycin is a monoprotic base (Table S1, pKa 9.5), so the fraction
    # unionised follows the single-pKa Henderson-Hasselbalch form; at pH 7.4
    # it is 0.0079, the factor Table S2 footnote e divides the QSAR
    # permeability by.
    fniFluid <- 1 / (1 + 10^(pka1 - phFluid))
    fniMass  <- 1 / (1 + 10^(pka1 - phMass))
    fniBlood <- 1 / (1 + 10^(pka1 - phBlood))

    # Gaohua 2015 Appendix S1 eq 1 with R = 8.314 Pa*m^3/mol/K.
    kaf <- henry / (8.314 * tempK)

    # Kanamycin is cleared renally only (Table S2: metabolic elimination 0).
    cl <- exp(lcl_renal)
    # =============== 3. Regional lung geometry ===============
    # PBR, mass and ELF split on the same proportions; a lobe and an airway
    # segment therefore each have one volume triple.
    vMassTotal   <- vLungTotal - vElfTotal - vPbrTotal
    vBloodLobe   <- vPbrTotal  * fSegLobe
    vMassLobe    <- vMassTotal * fSegLobe
    vFluidLobe   <- vElfTotal  * fSegLobe
    vBloodAirway <- vPbrTotal  * fSegAirway
    vMassAirway  <- vMassTotal * fSegAirway
    vFluidAirway <- vElfTotal  * fSegAirway

    vAirRt <- vAlvTotal * fAlvRt
    vAirRm <- vAlvTotal * fAlvRm
    vAirRl <- vAlvTotal * fAlvRl
    vAirLt <- vAlvTotal * fAlvLt
    vAirLl <- vAlvTotal * fAlvLl

    # Gaohua 2015 Appendix S1 eq 2: the fluid and air compartments are in immediate
    # equilibrium, so they are solved as one state with effective volume
    # V_AF = V_F + K_AF * V_A.
    vafRt <- vFluidLobe   + kaf * vAirRt
    vafRm <- vFluidLobe   + kaf * vAirRm
    vafRl <- vFluidLobe   + kaf * vAirRl
    vafLt <- vFluidLobe   + kaf * vAirLt
    vafLl <- vFluidLobe   + kaf * vAirLl
    vafLa <- vFluidAirway + kaf * vAirwayAir
    vafUa <- vFluidAirway + kaf * vAirwayAir

    # Permeation surface area (m^2) per segment.
    saRt <- saDeep * fAlvRt
    saRm <- saDeep * fAlvRm
    saRl <- saDeep * fAlvRl
    saLt <- saDeep * fAlvLt
    saLl <- saDeep * fAlvLl
    saLa <- saAirway / 2
    saUa <- saAirway / 2

    # Permeability-surface product, CL_PD = Peff * SA.  Peff is in cm/s and SA
    # in m^2, so the conversion to L/h is 1e4 cm^2/m^2 * 3600 s/h / 1000 mL/L
    # = 3.6e4.  Following the paper's single published "lung effective
    # permeability" per compound and its single in vitro-in vivo
    # extrapolation, the basal (mass-blood) product is the apical
    # (fluid-mass) one scaled by ratioPdBasal, which ships as 1.
    clpdRt <- peffLung * saRt * 3.6e4
    clpdRm <- peffLung * saRm * 3.6e4
    clpdRl <- peffLung * saRl * 3.6e4
    clpdLt <- peffLung * saLt * 3.6e4
    clpdLl <- peffLung * saLl * 3.6e4
    clpdLa <- peffLung * saLa * 3.6e4
    clpdUa <- peffLung * saUa * 3.6e4

    # Basal (mass-blood) permeability-surface products.
    clpdRtMb <- clpdRt * ratioPdBasal
    clpdRmMb <- clpdRm * ratioPdBasal
    clpdRlMb <- clpdRl * ratioPdBasal
    clpdLtMb <- clpdLt * ratioPdBasal
    clpdLlMb <- clpdLl * ratioPdBasal
    clpdLaMb <- clpdLa * ratioPdBasal
    clpdUaMb <- clpdUa * ratioPdBasal

    # Blood flow (L/h) and ventilation rate (L/h) per segment.
    qRt <- fqRt * qc
    qRm <- fqRm * qc
    qRl <- fqRl * qc
    qLt <- fqLt * qc
    qLl <- fqLl * qc
    qLa <- fqLa * qc
    qUa <- fqUa * qc

    ventTotal <- vqRatio * qc
    ventRt    <- fvRt * ventTotal
    ventRm    <- fvRm * ventTotal
    ventRl    <- fvRl * ventTotal
    ventLt    <- fvLt * ventTotal
    ventLl    <- fvLl * ventTotal

    # =============== 4. Compartment concentrations ===============
    Cc   <- central / vc
    Cvb  <- Cc * bp
    Cab  <- arterial / vArterial
    Cpbr <- lung_pbr / vPbrTotal

    Cblood_rt <- lung_rt_blood / vBloodLobe
    Cblood_rm <- lung_rm_blood / vBloodLobe
    Cblood_rl <- lung_rl_blood / vBloodLobe
    Cblood_lt <- lung_lt_blood / vBloodLobe
    Cblood_ll <- lung_ll_blood / vBloodLobe
    Cblood_la <- lung_la_blood / vBloodAirway
    Cblood_ua <- lung_ua_blood / vBloodAirway

    Cmass_rt <- lung_rt_mass / vMassLobe
    Cmass_rm <- lung_rm_mass / vMassLobe
    Cmass_rl <- lung_rl_mass / vMassLobe
    Cmass_lt <- lung_lt_mass / vMassLobe
    Cmass_ll <- lung_ll_mass / vMassLobe
    Cmass_la <- lung_la_mass / vMassAirway
    Cmass_ua <- lung_ua_mass / vMassAirway

    Celf_rt <- lung_rt_fluid / vafRt
    Celf_rm <- lung_rm_fluid / vafRm
    Celf_rl <- lung_rl_fluid / vafRl
    Celf_lt <- lung_lt_fluid / vafLt
    Celf_ll <- lung_ll_fluid / vafLl
    Celf_la <- lung_la_fluid / vafLa
    Celf_ua <- lung_ua_fluid / vafUa

    # Unbound, unionised driving concentrations.  Only unbound unionised drug
    # is passively permeable (Gaohua 2015 Methods "Passive permeability estimates"), so
    # every passive term is written on fu * (1/ionisation) * C.  fu in blood
    # follows from the definition of the blood:plasma ratio.
    fuBlood <- fup / bp

    uElf_rt   <- fniFluid * fuFluid * Celf_rt
    uElf_rm   <- fniFluid * fuFluid * Celf_rm
    uElf_rl   <- fniFluid * fuFluid * Celf_rl
    uElf_lt   <- fniFluid * fuFluid * Celf_lt
    uElf_ll   <- fniFluid * fuFluid * Celf_ll
    uElf_la   <- fniFluid * fuFluid * Celf_la
    uElf_ua   <- fniFluid * fuFluid * Celf_ua

    uMass_rt  <- fniMass * fuMass * Cmass_rt
    uMass_rm  <- fniMass * fuMass * Cmass_rm
    uMass_rl  <- fniMass * fuMass * Cmass_rl
    uMass_lt  <- fniMass * fuMass * Cmass_lt
    uMass_ll  <- fniMass * fuMass * Cmass_ll
    uMass_la  <- fniMass * fuMass * Cmass_la
    uMass_ua  <- fniMass * fuMass * Cmass_ua

    uBlood_rt <- fniBlood * fuBlood * Cblood_rt
    uBlood_rm <- fniBlood * fuBlood * Cblood_rm
    uBlood_rl <- fniBlood * fuBlood * Cblood_rl
    uBlood_lt <- fniBlood * fuBlood * Cblood_lt
    uBlood_ll <- fniBlood * fuBlood * Cblood_ll
    uBlood_la <- fniBlood * fuBlood * Cblood_la
    uBlood_ua <- fniBlood * fuBlood * Cblood_ua

    # =============== 5. ODE system ===============
    d/dt(depot) <- -ka * depot
    f(depot)    <- exp(lfdepot)

    # right lung top lobe (Gaohua 2015 Appendix S1 eq 3-5 / 6-8 written per segment).
    d/dt(lung_rt_fluid) <- ventRt * kaf * (Celf_la - Celf_rt) +
      clpdRt * (uMass_rt - uElf_rt) +
      clEffluxFm * fuMass * Cmass_rt - clUptakeFm * fuFluid * Celf_rt
    d/dt(lung_rt_mass) <- clpdRt * (uElf_rt - uMass_rt) +
      clUptakeFm * fuFluid * Celf_rt - clEffluxFm * fuMass * Cmass_rt +
      clpdRtMb * (uBlood_rt - uMass_rt) +
      clUptakeMb * fuBlood * Cblood_rt - clEffluxMb * fuMass * Cmass_rt -
      clMet * fuMass * Cmass_rt
    d/dt(lung_rt_blood) <- qRt * (Cpbr - Cblood_rt) +
      clpdRtMb * (uMass_rt - uBlood_rt) +
      clEffluxMb * fuMass * Cmass_rt - clUptakeMb * fuBlood * Cblood_rt

    # right lung middle lobe (Gaohua 2015 Appendix S1 eq 3-5 / 6-8 written per segment).
    d/dt(lung_rm_fluid) <- ventRm * kaf * (Celf_la - Celf_rm) +
      clpdRm * (uMass_rm - uElf_rm) +
      clEffluxFm * fuMass * Cmass_rm - clUptakeFm * fuFluid * Celf_rm
    d/dt(lung_rm_mass) <- clpdRm * (uElf_rm - uMass_rm) +
      clUptakeFm * fuFluid * Celf_rm - clEffluxFm * fuMass * Cmass_rm +
      clpdRmMb * (uBlood_rm - uMass_rm) +
      clUptakeMb * fuBlood * Cblood_rm - clEffluxMb * fuMass * Cmass_rm -
      clMet * fuMass * Cmass_rm
    d/dt(lung_rm_blood) <- qRm * (Cpbr - Cblood_rm) +
      clpdRmMb * (uMass_rm - uBlood_rm) +
      clEffluxMb * fuMass * Cmass_rm - clUptakeMb * fuBlood * Cblood_rm

    # right lung lower lobe (Gaohua 2015 Appendix S1 eq 3-5 / 6-8 written per segment).
    d/dt(lung_rl_fluid) <- ventRl * kaf * (Celf_la - Celf_rl) +
      clpdRl * (uMass_rl - uElf_rl) +
      clEffluxFm * fuMass * Cmass_rl - clUptakeFm * fuFluid * Celf_rl
    d/dt(lung_rl_mass) <- clpdRl * (uElf_rl - uMass_rl) +
      clUptakeFm * fuFluid * Celf_rl - clEffluxFm * fuMass * Cmass_rl +
      clpdRlMb * (uBlood_rl - uMass_rl) +
      clUptakeMb * fuBlood * Cblood_rl - clEffluxMb * fuMass * Cmass_rl -
      clMet * fuMass * Cmass_rl
    d/dt(lung_rl_blood) <- qRl * (Cpbr - Cblood_rl) +
      clpdRlMb * (uMass_rl - uBlood_rl) +
      clEffluxMb * fuMass * Cmass_rl - clUptakeMb * fuBlood * Cblood_rl

    # left lung top lobe (Gaohua 2015 Appendix S1 eq 3-5 / 6-8 written per segment).
    d/dt(lung_lt_fluid) <- ventLt * kaf * (Celf_la - Celf_lt) +
      clpdLt * (uMass_lt - uElf_lt) +
      clEffluxFm * fuMass * Cmass_lt - clUptakeFm * fuFluid * Celf_lt
    d/dt(lung_lt_mass) <- clpdLt * (uElf_lt - uMass_lt) +
      clUptakeFm * fuFluid * Celf_lt - clEffluxFm * fuMass * Cmass_lt +
      clpdLtMb * (uBlood_lt - uMass_lt) +
      clUptakeMb * fuBlood * Cblood_lt - clEffluxMb * fuMass * Cmass_lt -
      clMet * fuMass * Cmass_lt
    d/dt(lung_lt_blood) <- qLt * (Cpbr - Cblood_lt) +
      clpdLtMb * (uMass_lt - uBlood_lt) +
      clEffluxMb * fuMass * Cmass_lt - clUptakeMb * fuBlood * Cblood_lt

    # left lung lower lobe (Gaohua 2015 Appendix S1 eq 3-5 / 6-8 written per segment).
    d/dt(lung_ll_fluid) <- ventLl * kaf * (Celf_la - Celf_ll) +
      clpdLl * (uMass_ll - uElf_ll) +
      clEffluxFm * fuMass * Cmass_ll - clUptakeFm * fuFluid * Celf_ll
    d/dt(lung_ll_mass) <- clpdLl * (uElf_ll - uMass_ll) +
      clUptakeFm * fuFluid * Celf_ll - clEffluxFm * fuMass * Cmass_ll +
      clpdLlMb * (uBlood_ll - uMass_ll) +
      clUptakeMb * fuBlood * Cblood_ll - clEffluxMb * fuMass * Cmass_ll -
      clMet * fuMass * Cmass_ll
    d/dt(lung_ll_blood) <- qLl * (Cpbr - Cblood_ll) +
      clpdLlMb * (uMass_ll - uBlood_ll) +
      clEffluxMb * fuMass * Cmass_ll - clUptakeMb * fuBlood * Cblood_ll

    # Lower airway.  The LA is the ventilation hub: it exchanges air with the
    # UA and with each of the five lobes (Gaohua 2015 Appendix S1 eq 21).
    d/dt(lung_la_fluid) <- ventTotal * kaf * (Celf_ua - Celf_la) +
      ventRt * kaf * (Celf_rt - Celf_la) +
      ventRm * kaf * (Celf_rm - Celf_la) +
      ventRl * kaf * (Celf_rl - Celf_la) +
      ventLt * kaf * (Celf_lt - Celf_la) +
      ventLl * kaf * (Celf_ll - Celf_la) +
      clpdLa * (uMass_la - uElf_la) +
      clEffluxFm * fuMass * Cmass_la - clUptakeFm * fuFluid * Celf_la
    d/dt(lung_la_mass) <- clpdLa * (uElf_la - uMass_la) +
      clUptakeFm * fuFluid * Celf_la - clEffluxFm * fuMass * Cmass_la +
      clpdLaMb * (uBlood_la - uMass_la) +
      clUptakeMb * fuBlood * Cblood_la - clEffluxMb * fuMass * Cmass_la -
      clMet * fuMass * Cmass_la
    d/dt(lung_la_blood) <- qLa * (Cpbr - Cblood_la) +
      clpdLaMb * (uMass_la - uBlood_la) +
      clEffluxMb * fuMass * Cmass_la - clUptakeMb * fuBlood * Cblood_la

    # Upper airway.  Inhaled air carries no drug and exhaled air leaves at the
    # UA air concentration (Gaohua 2015 Appendix S1 eq 24); the UA is perfused directly by
    # arterial blood and drains to the venous side (eq 26 and eq 8).
    d/dt(lung_ua_fluid) <- ventTotal * kaf * (0 - Celf_ua) +
      ventTotal * kaf * (Celf_la - Celf_ua) +
      clpdUa * (uMass_ua - uElf_ua) +
      clEffluxFm * fuMass * Cmass_ua - clUptakeFm * fuFluid * Celf_ua
    d/dt(lung_ua_mass) <- clpdUa * (uElf_ua - uMass_ua) +
      clUptakeFm * fuFluid * Celf_ua - clEffluxFm * fuMass * Cmass_ua +
      clpdUaMb * (uBlood_ua - uMass_ua) +
      clUptakeMb * fuBlood * Cblood_ua - clEffluxMb * fuMass * Cmass_ua -
      clMet * fuMass * Cmass_ua
    d/dt(lung_ua_blood) <- qUa * (Cab - Cblood_ua) +
      clpdUaMb * (uMass_ua - uBlood_ua) +
      clEffluxMb * fuMass * Cmass_ua - clUptakeMb * fuBlood * Cblood_ua

    # Pulmonary blood reservoir (Gaohua 2015 Appendix S1 eq 6).  Venous blood enters at
    # cardiac output; the five lobes and the lower airway return to it.  (The
    # recovered equation prints C_RLB in the RM, RT and LT terms, which is a
    # transcription artefact of the embedded-equation decode: each segment
    # returns its own blood concentration.)
    d/dt(lung_pbr) <- qc * (Cvb - Cpbr) +
      qRt * (Cblood_rt - Cpbr) +
      qRm * (Cblood_rm - Cpbr) +
      qRl * (Cblood_rl - Cpbr) +
      qLt * (Cblood_lt - Cpbr) +
      qLl * (Cblood_ll - Cpbr)  +
      qLa * (Cblood_la - Cpbr)

    # Arterial blood (Gaohua 2015 Appendix S1 eq 7).
    d/dt(arterial) <- qc * (Cpbr - Cab)

    # Systemic compartment (Gaohua 2015 Appendix S1 eq 8, reduced).  In the published
    # model this is venous blood fed by twelve perfusion-limited tissues; here
    # the twelve tissues and venous blood are lumped into one well-stirred
    # compartment of plasma-referenced volume vc, so arterial blood perfusing
    # the body returns at the systemic blood concentration.  Elimination is
    # the compound file's total systemic clearance on plasma concentration.
    d/dt(central) <- ka * depot +
      (qc - qUa) * Cab + qUa * Cblood_ua - qc * Cvb -
      cl * Cc

    # =============== 6. Observation ===============
    # Deterministic forward-simulation model: the paper reports no IIV and no
    # residual-error model.  Cc is the systemic plasma concentration; the
    # per-segment ELF (Celf_*) and tissue-mass (Cmass_*) concentrations above
    # are the model's reported outputs.
    Cc ~ add(addSd)
  })
}
