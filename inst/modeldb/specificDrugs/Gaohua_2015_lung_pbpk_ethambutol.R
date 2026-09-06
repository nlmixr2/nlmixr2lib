Gaohua_2015_lung_pbpk_ethambutol <- function() {
  description <- paste(
    "PBPK (permeability-limited lung, 25 ODEs; Simcyp V14 R1). Ethambutol",
    "plasma, epithelial lining fluid and lung tissue concentrations in the",
    "multicompartment permeability-limited lung PBPK model of Gaohua et al.",
    "(Simcyp V14 R1).  The lung is divided into seven segments - upper and",
    "lower airways plus the five lobes - each carrying four compartments:",
    "pulmonary capillary blood, tissue mass, fluid (mucus and epithelial",
    "lining fluid) and alveolar air, giving the 21 lung ODEs of Appendix S1",
    "plus the pulmonary blood reservoir and arterial blood.  Passive",
    "permeation of unbound unionised drug crosses the apical (fluid-mass)",
    "and basal (mass-blood) membranes, with Henderson-Hasselbalch ionisation",
    "at the fluid (pH 6.6), mass (pH 6.69) and blood (pH 7.4) pH values;",
    "uptake and efflux transporters and lung metabolism are exposed and",
    "fixed at zero, which is the paper's base case.  DEVIATION FROM THE",
    "PUBLISHED ARCHITECTURE: the paper embeds this lung model in a Simcyp",
    "full PBPK whose twelve perfusion-limited tissues need per-tissue Kp",
    "values that are not published for any compound, so the systemic side is",
    "reduced here to one well-stirred compartment at the compound file's own",
    "Vss, ka, fa and clearance.  Deterministic typical-value simulation",
    "model: the paper reports no IIV and no residual-error model.  See the",
    "validation vignette for the reduction's accuracy against the paper's own",
    "predicted plasma profiles."
  )
  reference <- paste(
    "Gaohua L, Wedagedera J, Small BG, Almond L, Romero K, Hermann D, Hanna",
    "D, Jamei M, Gardner I. Development of a Multicompartment",
    "Permeability-Limited Lung PBPK Model and Its Application in Predicting",
    "Pulmonary Pharmacokinetics of Antituberculosis Drugs. CPT Pharmacometrics",
    "Syst Pharmacol. 2015;4(10):605-613. doi:10.1002/psp4.12034.",
    "Reference physiology (cardiac output 356 L/h, lung volume 0.53 L, blood",
    "volume 5.75 L) is taken from the upstream full-body PBPK the Appendix",
    "defers to: Jamei M, Bajot F, Neuhoff S, Barter Z, Yang J,",
    "Rostami-Hodjegan A, Rowland-Yeo K. Clin Pharmacokinet. 2014;53:73-87,",
    "Electronic Supplementary Material 1 Tables S1-S3.",
    sep = " "
  )
  vignette <- "Gaohua_2015_lung_pbpk"
  units    <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # The lung segments and the pulmonary blood reservoir are anatomical states
  # of this paper's lung model and do not generalise to other extractions;
  # `arterial` is the systemic arterial blood pool of Appendix S1 eq 7.
  paper_specific_compartments        <- c("arterial")
  paper_specific_compartment_pattern <- "^lung_"

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Scales the systemic distribution volume only (Vss is reported in",
        "L/kg in Supplementary Table S2).  Clearance and the lung physiology",
        "are absolute values for the reference adult and are not weight",
        "scaled; the paper reports no allometric relationship.  Use 70 kg to",
        "reproduce the simulations in this package's vignette.",
        sep = " "
      ),
      source_name        = "WT"
    )
  )

  compartmentData <- list(
    depot        = list(analyte = "ethambutol", units = "mg", specimen = "administration site", verified = TRUE),
    central      = list(analyte = "ethambutol", units = "mg", specimen = "plasma", verified = TRUE),
    arterial     = list(analyte = "ethambutol", units = "mg", specimen = "whole blood", verified = TRUE),
    lung_pbr     = list(analyte = "ethambutol", units = "mg", specimen = "whole blood", verified = TRUE),
    lung_rt_fluid = list(analyte = "ethambutol", units = "mg", specimen = "epithelial lining fluid", verified = TRUE),
    lung_rt_mass  = list(analyte = "ethambutol", units = "mg", specimen = "tissue", verified = TRUE),
    lung_rt_blood = list(analyte = "ethambutol", units = "mg", specimen = "whole blood", verified = TRUE),
    lung_rm_fluid = list(analyte = "ethambutol", units = "mg", specimen = "epithelial lining fluid", verified = TRUE),
    lung_rm_mass  = list(analyte = "ethambutol", units = "mg", specimen = "tissue", verified = TRUE),
    lung_rm_blood = list(analyte = "ethambutol", units = "mg", specimen = "whole blood", verified = TRUE),
    lung_rl_fluid = list(analyte = "ethambutol", units = "mg", specimen = "epithelial lining fluid", verified = TRUE),
    lung_rl_mass  = list(analyte = "ethambutol", units = "mg", specimen = "tissue", verified = TRUE),
    lung_rl_blood = list(analyte = "ethambutol", units = "mg", specimen = "whole blood", verified = TRUE),
    lung_lt_fluid = list(analyte = "ethambutol", units = "mg", specimen = "epithelial lining fluid", verified = TRUE),
    lung_lt_mass  = list(analyte = "ethambutol", units = "mg", specimen = "tissue", verified = TRUE),
    lung_lt_blood = list(analyte = "ethambutol", units = "mg", specimen = "whole blood", verified = TRUE),
    lung_ll_fluid = list(analyte = "ethambutol", units = "mg", specimen = "epithelial lining fluid", verified = TRUE),
    lung_ll_mass  = list(analyte = "ethambutol", units = "mg", specimen = "tissue", verified = TRUE),
    lung_ll_blood = list(analyte = "ethambutol", units = "mg", specimen = "whole blood", verified = TRUE),
    lung_la_fluid = list(analyte = "ethambutol", units = "mg", specimen = "epithelial lining fluid", verified = TRUE),
    lung_la_mass  = list(analyte = "ethambutol", units = "mg", specimen = "tissue", verified = TRUE),
    lung_la_blood = list(analyte = "ethambutol", units = "mg", specimen = "whole blood", verified = TRUE),
    lung_ua_fluid = list(analyte = "ethambutol", units = "mg", specimen = "epithelial lining fluid", verified = TRUE),
    lung_ua_mass  = list(analyte = "ethambutol", units = "mg", specimen = "tissue", verified = TRUE),
    lung_ua_blood = list(analyte = "ethambutol", units = "mg", specimen = "whole blood", verified = TRUE)
  )

  population <- list(
    species     = "human",
    n_subjects  = 540L,
    n_studies   = 2L,
    age_range   = "25-60 years",
    sex_female_pct = 42.9,
    disease_state  = "virtual healthy North European Caucasian adults",
    dose_range     = "25 mg/kg single oral dose (plasma verification); 15 mg/kg once daily for 5 days (lung)",
    regions        = "North European Caucasian (Simcyp virtual population)",
    notes       = paste(
      "Supplementary Appendix S2: the plasma-verification simulation used 140",
      "North European Caucasian subjects (age 25-60, 42.9% female) dosed 25",
      "mg/kg and followed for 72 h, compared with observed data from Peloquin",
      "et al. 1999.  The lung simulation used 400 subjects (age 25-50, 50%",
      "female) dosed 15 mg/kg once daily for 5 days, compared with Conte et",
      "al. 2001 (Supplementary Table S1).  These are virtual subjects",
      "generated by the Simcyp population library, not enrolled participants.",
      sep = " "
    )
  )

  ini({
    # ================= System: reference physiology =================
    # Cardiac output and the total lung / blood volumes come from the upstream
    # Simcyp full-body PBPK the Appendix defers to (Jamei 2014 ESM1); every
    # regional split below is printed in the main text under "Model
    # parameters".
    qc         <- fixed(356);   label("Cardiac output (L/h)")                              # Jamei 2014 ESM1 Table S3 ("On average the cardiac output (CO) is 356 (L/h)")
    vLungTotal <- fixed(0.53);  label("Total lung volume (L)")                             # Jamei 2014 ESM1 Table S1 (Lungs 0.53 L)
    vArterial  <- fixed(1.955); label("Arterial blood volume (L)")                         # Jamei 2014 ESM1 Table S1 (blood 5.75 L; "venous blood volume is assumed to be 66% of the total blood volume") -> 0.34 * 5.75
    vPbrTotal  <- fixed(0.089); label("Total pulmonary blood reservoir volume (L)")        # Methods "Tissue volumes" (mean 89 mL, CV 20%, literature range 63-114 mL)
    vElfTotal  <- fixed(0.025); label("Total epithelial lining fluid volume (L)")          # Methods "Tissue volumes" (mean 25 mL, CV 20%, literature range 10-40 mL)
    vAlvTotal  <- fixed(5.6);   label("Total alveolar air volume (L)")                     # Methods "Distribution of alveolar volume" (5.6 L, CV 10%)
    vAirwayAir <- fixed(0.05);  label("Air volume of one airway segment (L)")              # Methods "Upper/lower airways" (UA and LA both 0.05 L)

    # PBR, lung mass and ELF are all split on the same proportions: 16.67% to
    # each of the five lobes, the remaining 16.67% shared equally by UA and LA.
    fSegLobe   <- fixed(0.1666667);  label("Fraction of PBR / mass / ELF in each lung lobe")     # Methods "Tissue volumes" (each of the five lobes receives 16.67% of the total PBR)
    fSegAirway <- fixed(0.08333333); label("Fraction of PBR / mass / ELF in each airway segment")# Methods "Tissue volumes" (remaining 16.67% split equally between the UA and LA segments)

    # Alveolar air volume by lobe (Methods "Distribution of alveolar volume").
    fAlvRt <- fixed(0.1578947); label("Fraction of alveolar air volume, right top lobe")    # Methods (3/19 RT)
    fAlvRm <- fixed(0.1052632); label("Fraction of alveolar air volume, right middle lobe") # Methods (2/19 RM)
    fAlvRl <- fixed(0.2631579); label("Fraction of alveolar air volume, right lower lobe")  # Methods (5/19 RL)
    fAlvLt <- fixed(0.2105263); label("Fraction of alveolar air volume, left top lobe")     # Methods (4/19 LT)
    fAlvLl <- fixed(0.2631579); label("Fraction of alveolar air volume, left lower lobe")   # Methods (5/19 LL)

    # Pulmonary blood flow as a fraction of cardiac output (Methods "Blood
    # flow rate"); the upper airway is perfused from arterial blood instead.
    fqRt <- fixed(0.086); label("Pulmonary blood flow fraction, right top lobe")    # Methods "Blood flow rate" (8.6% RT)
    fqRm <- fixed(0.118); label("Pulmonary blood flow fraction, right middle lobe") # Methods "Blood flow rate" (11.8% RM)
    fqRl <- fixed(0.310); label("Pulmonary blood flow fraction, right lower lobe")  # Methods "Blood flow rate" (31.0% RL)
    fqLt <- fixed(0.086); label("Pulmonary blood flow fraction, left top lobe")     # Methods "Blood flow rate" (8.6% LT)
    fqLl <- fixed(0.349); label("Pulmonary blood flow fraction, left lower lobe")   # Methods "Blood flow rate" (34.9% LL)
    fqLa <- fixed(0.050); label("Pulmonary blood flow fraction, lower airway")      # Methods "Blood flow rate" (5.0% LA)
    fqUa <- fixed(0.025); label("Arterial blood flow fraction, upper airway")       # Methods "Blood flow rate" (2.5% of arterial blood circulates in the UA)

    # Alveolar ventilation.  The total is set from the ventilation/perfusion
    # ratio and then partitioned by lobe (Methods "Ventilation/perfusion
    # distribution").
    vqRatio <- fixed(1.0);   label("Ventilation/perfusion ratio (unitless)")   # Methods "Ventilation/perfusion distribution" (geometric mean = 1.0)
    fvRt    <- fixed(0.149); label("Ventilation fraction, right top lobe")     # Methods (14.9% RT)
    fvRm    <- fixed(0.124); label("Ventilation fraction, right middle lobe")  # Methods (12.4% RM)
    fvRl    <- fixed(0.259); label("Ventilation fraction, right lower lobe")   # Methods (25.9% RL)
    fvLt    <- fixed(0.148); label("Ventilation fraction, left top lobe")      # Methods (14.8% LT)
    fvLl    <- fixed(0.320); label("Ventilation fraction, left lower lobe")    # Methods (32.0% LL)

    # Absorption (permeation) surface area: the deep lung is apportioned by
    # relative alveolar volume, the airways are split 50:50 (Methods
    # "Absorption area").
    saDeep   <- fixed(140); label("Deep-lung permeation surface area (m^2)")   # Methods "Absorption area" (~140 m2, CV 30%, assigned to each lobe by relative alveolar volume)
    saAirway <- fixed(1.5); label("Upper + lower airway surface area (m^2)")   # Methods "Absorption area" (1.5 m2, CV 50%, for the sum of the UA and LA, split 50:50)

    # Compartment pH used in the Henderson-Hasselbalch ionisation terms.
    phFluid <- fixed(6.60); label("pH of the lung fluid (ELF) compartment")    # Methods "pH" (pH about 6.6 in healthy individuals)
    phMass  <- fixed(6.69); label("pH of the lung tissue mass compartment")    # Methods "pH" (pH in the lung mass is around 6.69 +/- 0.07)
    phBlood <- fixed(7.40); label("pH of blood / plasma")                      # Methods "pH" (arterial blood pH between 7.38 and 7.43); the ionisation corrections in Results are quoted at pH 7.4

    # Lung metabolism and the four transporter clearances are the paper's base
    # case: passive diffusion only.  They are exposed so a user can reproduce
    # the published sensitivity analyses.
    clMet      <- fixed(0); label("Metabolic clearance in the lung mass compartment (L/h)")            # Methods "Parameterization of the multicompartment lung model" (for all compounds, lung metabolism was assumed to be negligible, CLmet = 0)
    clUptakeFm <- fixed(0); label("Apical (fluid->mass) uptake transporter clearance (L/h)")           # Methods (initial simulations considered distribution by passive diffusion only)
    clEffluxFm <- fixed(0); label("Apical (mass->fluid) efflux transporter clearance (L/h)")           # Methods; Figure 5 varies this over 0, 0.06, 0.6, 6 and 60 L/h
    clUptakeMb <- fixed(0); label("Basal (blood->mass) uptake transporter clearance (L/h)")            # Methods; Supplementary Figure 4 varies this over 0, 0.6 and 60 L/h
    clEffluxMb <- fixed(0); label("Basal (mass->blood) efflux transporter clearance (L/h)")            # Methods (initial simulations considered distribution by passive diffusion only)

    # Appendix S1 requires TWO permeability-surface products per segment --
    # apical (fluid-mass, CL_PD,FM) and basal (mass-blood, CL_PD,MB) -- but
    # Table S2 publishes ONE "lung effective permeability" per compound and
    # the Methods extrapolate it using only the ABSORPTION surface areas.  The
    # paper never states whether the two are equal.  The shipped default of 1
    # makes them equal, which is the reading most consistent with a single
    # published permeability and a single in vitro-in vivo extrapolation.  It
    # is exposed as a parameter because it is an unpublished structural
    # choice, and because it is the one assumption that measurably misses the
    # paper's own ethambutol result (see the validation vignette).
    ratioPdBasal <- fixed(1); label("Basal:apical permeability-surface ratio (unitless)")               # NOT PUBLISHED: Appendix S1 needs CL_PD,MB and CL_PD,FM separately; Table S2 gives one permeability

    # Air:fluid partition coefficient, K_AF = H / (R * T) (Appendix S1 eq 1).
    # Henry's constant H was "predicted using a QSAR approach" and is not
    # printed for any compound.  All three compounds here are non-volatile
    # crystalline solids, for which H -> 0 and therefore K_AF -> 0: the
    # alveolar-air compartment holds no drug and the ventilation terms of
    # equations 3, 21 and 24 drop out.  Exposed as a parameter so the
    # ventilation limb can be switched on for a volatile compound.
    kaf <- fixed(0); label("Air:fluid partition coefficient (unitless)")  # Appendix S1 eq 1; H not printed, taken to the non-volatile limit H -> 0

    # ================= Compound file: ethambutol =================
    # Supplementary Table S2 unless noted.
    bp       <- fixed(1.3);   label("Blood:plasma concentration ratio")               # Table S2 "B/P"
    fup      <- fixed(0.75);  label("Unbound fraction in plasma")                     # Table S2 "fu p" (main plasma binding protein albumin)
    pka1     <- fixed(9.55);  label("First pKa (diprotic base)")                      # Table S2 "pKa" (9.55, 6.5); compound type "Diprotic Base"
    pka2     <- fixed(6.5);   label("Second pKa (diprotic base)")                     # Table S2 "pKa" (9.55, 6.5)
    lka      <- fixed(log(0.25)); label("First-order absorption rate constant (1/h)") # Table S2 "ka" (0.25 1/h)
    lfdepot  <- fixed(log(0.57)); label("Fraction of the dose absorbed (unitless)")   # Table S2 "fa" (0.57)
    lvc      <- fixed(log(1.23)); label("Systemic distribution volume (L/kg)")        # Table S2 "V ss" (1.23 L/kg)
    lcl_renal  <- fixed(log(25.55)); label("Renal clearance (L/h)")                   # Table S2 "CL R" (25.55 L/h)
    lcl_nonren <- fixed(log(4.35));  label("Additional systemic clearance (L/h)")     # Table S2 "Additional systemic CL" (4.35 L/h)

    peffLung <- fixed(0.479e-4); label("Lung effective permeability (cm/s)")            # Table S2 "Lung effective permeability" (0.479 x 10^-4 cm/s); Results 'Data sources': QSAR-predicted 3.0 x 10^-7 cm/s divided by fni = 0.0062 at pH 7.4
    fuMass   <- fixed(0.451); label("Unbound fraction in the lung tissue mass")            # Table S2 "Lung fu mass"; Methods reports the same value for ethambutol
    fuFluid  <- fixed(1); label("Unbound fraction in the lung fluid (ELF)")           # Table S2 "fu fluid"

    # No residual-error model is reported: the paper is a deterministic
    # forward simulation.  An additive term fixed at zero keeps the model
    # solvable by nlmixr2 without inventing a variance.
    addSd <- fixed(0); label("Additive residual error (mg/L)")                        # not reported; the paper presents predictions, not a fit
  })

  model({
    # =============== 1. Reference-individual scaling ===============
    # The systemic reduction (see description): one well-stirred systemic
    # compartment carrying the whole body except the lung, at the compound
    # file's own Vss, with the lung layer solved exactly as published.
    vc <- exp(lvc) * WT
    ka <- exp(lka)

    # =============== 2. Ionisation and total clearance ===============
    # Ethambutol is a diprotic base (Table S2, pKa 9.55 and 6.5).  At pH 7.4
    # this form gives fni = 0.00625, reproducing the paper's stated 0.0062,
    # and Table S2's lung effective permeability 0.479e-4 cm/s is recovered as
    # the QSAR 3.0e-7 cm/s divided by that same 0.0062.
    fniFluid <- 1 / (1 + 10^(pka1 - phFluid) + 10^(pka1 + pka2 - 2 * phFluid))
    fniMass  <- 1 / (1 + 10^(pka1 - phMass)  + 10^(pka1 + pka2 - 2 * phMass))
    fniBlood <- 1 / (1 + 10^(pka1 - phBlood) + 10^(pka1 + pka2 - 2 * phBlood))

    # Total systemic clearance is the sum of the two clearance terms in the
    # compound file.
    cl <- exp(lcl_renal) + exp(lcl_nonren)

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

    # Appendix S1 eq 2: the fluid and air compartments are in immediate
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
    # is passively permeable (Methods "Passive permeability estimates"), so
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

    # right lung top lobe (Appendix S1 eq 3-5 / 6-8 written per segment).
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

    # right lung middle lobe (Appendix S1 eq 3-5 / 6-8 written per segment).
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

    # right lung lower lobe (Appendix S1 eq 3-5 / 6-8 written per segment).
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

    # left lung top lobe (Appendix S1 eq 3-5 / 6-8 written per segment).
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

    # left lung lower lobe (Appendix S1 eq 3-5 / 6-8 written per segment).
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
    # UA and with each of the five lobes (Appendix S1 eq 21).
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
    # UA air concentration (Appendix S1 eq 24); the UA is perfused directly by
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

    # Pulmonary blood reservoir (Appendix S1 eq 6).  Venous blood enters at
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

    # Arterial blood (Appendix S1 eq 7).
    d/dt(arterial) <- qc * (Cpbr - Cab)

    # Systemic compartment (Appendix S1 eq 8, reduced).  In the published
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
