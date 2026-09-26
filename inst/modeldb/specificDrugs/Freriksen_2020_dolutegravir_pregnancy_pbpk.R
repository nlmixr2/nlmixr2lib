Freriksen_2020_dolutegravir_pregnancy_pbpk <- function() {
  description <- paste(
    "PBPK (whole-body pregnancy p-PBPK with a fetoplacental unit, hand-written",
    "in Berkeley Madonna). Oral dolutegravir in a typical pregnant woman in the",
    "third trimester (gestational age 34 weeks), with maternal and fetal",
    "plasma and amniotic-fluid concentrations. The 15 maternal ODE states are",
    "the oral depot (gut lumen), gut wall, lung, adipose, bone, brain, heart,",
    "kidney, muscle, skin, spleen, liver, rest of body, arterial and venous",
    "blood, all perfusion-rate limited with Simcyp-predicted Kp values; the",
    "liver is well-stirred with an IVIVE UGT1A1 + CYP3A4 intrinsic clearance",
    "multiplied by an empirical scaling factor of 65 fitted to the Min 2010",
    "healthy-volunteer data. The placenta is a barrier: bidirectional",
    "clearances measured in ex vivo dual-side cotyledon perfusion, corrected",
    "for perfusion-buffer protein binding and scaled by 30 cotyledons, link",
    "maternal arterial blood to a 3-state fetal model (fetal blood, rest of",
    "fetal body, amniotic fluid). Maternal and fetal physiology follow",
    "gestational-age regressions (covariate GA). Deterministic typical-value",
    "forward simulation: the paper has no IIV and no residual-error model."
  )
  reference <- paste(
    "Freriksen JJM, Schalkwijk S, Colbers AP, Abduljalil K, Russel FGM,",
    "Burger DM, Greupink R. Assessment of Maternal and Fetal Dolutegravir",
    "Exposure by Integrating Ex Vivo Placental Perfusion Data and",
    "Physiologically-Based Pharmacokinetic Modeling.",
    "Clin Pharmacol Ther. 2020;107(6):1352-1361. doi:10.1002/cpt.1748.",
    "Model code: Supplementary File S2 (Berkeley Madonna script); parameter",
    "tables S1-S2."
  )
  vignette <- "Freriksen_2020_dolutegravir_pregnancy_pbpk"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  paper_specific_compartments <- c("rest", "blood_fet", "rest_fet", "amniotic")

  compartmentData <- list(
    depot = list(analyte = "dolutegravir", units = "mg", specimen = "administration site", verified = TRUE),
    gut = list(analyte = "dolutegravir", units = "mg", specimen = "tissue", verified = TRUE),
    lung = list(analyte = "dolutegravir", units = "mg", specimen = "tissue", verified = TRUE),
    adipose = list(analyte = "dolutegravir", units = "mg", specimen = "tissue", verified = TRUE),
    bone = list(analyte = "dolutegravir", units = "mg", specimen = "tissue", verified = TRUE),
    brain = list(analyte = "dolutegravir", units = "mg", specimen = "tissue", verified = TRUE),
    heart = list(analyte = "dolutegravir", units = "mg", specimen = "tissue", verified = TRUE),
    kidney = list(analyte = "dolutegravir", units = "mg", specimen = "tissue", verified = TRUE),
    muscle = list(analyte = "dolutegravir", units = "mg", specimen = "tissue", verified = TRUE),
    skin = list(analyte = "dolutegravir", units = "mg", specimen = "tissue", verified = TRUE),
    spleen = list(analyte = "dolutegravir", units = "mg", specimen = "tissue", verified = TRUE),
    liver = list(analyte = "dolutegravir", units = "mg", specimen = "tissue", verified = TRUE),
    rest = list(analyte = "dolutegravir", units = "mg", specimen = "tissue", verified = TRUE),
    arterial = list(analyte = "dolutegravir", units = "mg", specimen = "whole blood", verified = TRUE),
    venous = list(analyte = "dolutegravir", units = "mg", specimen = "whole blood", verified = TRUE),
    blood_fet = list(analyte = "dolutegravir", units = "mg", specimen = "whole blood", verified = TRUE),
    rest_fet = list(analyte = "dolutegravir", units = "mg", specimen = "tissue", verified = TRUE),
    amniotic = list(analyte = "dolutegravir", units = "mg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list(
    GA = list(
      description = "Gestational age",
      units = "weeks",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Drives every maternal and fetal physiology regression of Table S2",
        "(body weight, cardiac output, CYP3A4 / UGT1A1 induction, fetal",
        "volume, fetal blood volume, fetal cardiac output, fraction of fetal",
        "cardiac output to the placenta, amniotic-fluid volume). Held",
        "constant over a simulation, as in the deposited script. The paper",
        "simulated and validated GA = 34 weeks only, and the protein-binding",
        "inputs fu and fu_fet are the 34-week values (they are not GA",
        "regressions), so other GA values extrapolate. The fetal-blood",
        "regression is negative below ~18.7 weeks: third trimester only."
      ),
      source_name = "GA"
    )
  )

  population <- list(
    species = "human",
    n_subjects = NA_integer_,
    n_studies = 2L,
    age_range = "not reported (typical adult pregnant woman)",
    weight_range = "73.7 kg at 34 weeks of gestation (Table S2 regression)",
    sex_female_pct = 100,
    disease_state = paste(
      "Pregnant women living with HIV-1 in the third trimester, and their",
      "fetuses. Maternal predictions were compared with mean third-trimester",
      "profiles from the PANNA network (n = 15) and IMPAACT P1026s (n = 28);",
      "fetal predictions with cord-to-maternal plasma ratios (PANNA n = 10,",
      "Rimawi n = 3, IMPAACT n = 23) and one amniotic-fluid measurement."
    ),
    dose_range = "Dolutegravir 50 mg orally once daily (7 doses); healthy-volunteer calibration used single doses of 2-100 mg (Min 2010)",
    regions = "Netherlands / Europe (PANNA), United States (IMPAACT)",
    notes = paste(
      "A single typical individual: covariates were matched to the mean",
      "values of the clinical studies. The placental transfer inputs come",
      "from n = 3 ex vivo dual-side cotyledon perfusions in each direction",
      "using healthy term placentas. The underlying nonpregnant",
      "healthy-volunteer PBPK (fu = 0.0041, nonpregnant physiology) was used",
      "only to fit the liver scaling factor and is not separately specified."
    )
  )

  ini({
    # Every value is a fixed input of the deposited Berkeley Madonna script
    # (Supplementary File S2); nothing was estimated except SF and A_SF,
    # which the authors fitted by eye against clinical data and then fixed.

    # Absorption (File S2 'ABSORPTION'; Table S1)
    lka <- fixed(log(2.25)); label("First-order absorption rate constant (1/h)") # Table S1 'Ka 2.25 h-1'; File S2 Ka = 2.25
    lfdepot <- fixed(log(0.804)); label("Fraction available from dosage form, Fabs (unitless)") # Table S1 'Fabs 0.804'; File S2 Fabs = 0.804

    # Maternal tissue:plasma partition coefficients (Simcyp v17, Kp scalar = 1)
    lkp_lung <- fixed(log(0.2237)); label("Lung tissue:plasma partition coefficient (unitless)") # Table S1 Kplu 0.2237
    lkp_adipose <- fixed(log(0.1163)); label("Adipose tissue:plasma partition coefficient (unitless)") # Table S1 Kpad 0.1163
    lkp_bone <- fixed(log(0.1953)); label("Bone tissue:plasma partition coefficient (unitless)") # Table S1 Kpbo 0.1953
    lkp_brain <- fixed(log(0.1395)); label("Brain tissue:plasma partition coefficient (unitless)") # Table S1 Kpbr 0.1395
    lkp_heart <- fixed(log(0.1825)); label("Heart tissue:plasma partition coefficient (unitless)") # Table S1 Kphe 0.1825
    lkp_kidney <- fixed(log(0.1674)); label("Kidney tissue:plasma partition coefficient (unitless)") # Table S1 Kpki 0.1674
    lkp_muscle <- fixed(log(0.0730)); label("Muscle tissue:plasma partition coefficient (unitless)") # Table S1 Kpmu 0.0730
    lkp_skin <- fixed(log(0.3195)); label("Skin tissue:plasma partition coefficient (unitless)") # Table S1 Kpsk 0.3195
    lkp_spleen <- fixed(log(0.1355)); label("Spleen tissue:plasma partition coefficient (unitless)") # Table S1 Kpspl 0.1355
    lkp_gut <- fixed(log(0.2294)); label("Gut wall tissue:plasma partition coefficient (unitless)") # Table S1 Kpgut 0.2294
    lkp_liver <- fixed(log(0.1447)); label("Liver tissue:plasma partition coefficient (unitless)") # Table S1 Kpliv 0.1447
    lkp_rest <- fixed(log(0.1752)); label("Rest-of-body tissue:plasma partition coefficient (unitless)") # Table S1 Kpre 0.1752 (average of all tissues)
    lkp_rest_fet <- fixed(log(0.1752)); label("Rest-of-fetal-body tissue:plasma partition coefficient (unitless)") # Table S1 KprbF 0.1752 (average of all tissues)

    # Binding (Table S1; pregnant-woman and fetal values at 34 weeks)
    fu <- fixed(0.0049); label("Maternal fraction unbound in plasma at 34 weeks (unitless)") # Table S1 'Pregnant women 0.0049'; File S2 fup_m
    fu_fet <- fixed(0.0036); label("Fetal fraction unbound in plasma (unitless)") # Table S1 'Fetus 0.0036'; File S2 fup_f
    bpr <- fixed(0.535); label("Maternal blood:plasma concentration ratio (unitless)") # Table S1 'Pregnant women 0.535'; File S2 BP
    bpr_fet <- fixed(0.441); label("Fetal blood:plasma concentration ratio (unitless)") # Table S1 'Fetus 0.441'; File S2 BP_f

    # Hepatic clearance (IVIVE from recombinant enzymes, Reese 2013)
    clint_ugt1a1 <- fixed(2.7); label("Apparent in vitro UGT1A1 intrinsic clearance (uL/min/mg)") # Table S1 CLint_met_appUGT1A1 2.7
    clint_cyp3a4 <- fixed(0.56); label("Apparent in vitro CYP3A4 intrinsic clearance (uL/min/mg)") # Table S1 CLint_met_appCYP3A4 0.56
    fumic <- fixed(0.817); label("Fraction unbound in the in vitro incubation (unitless)") # Table S1 fumic 0.817
    mppgl <- fixed(39.79067); label("Microsomal protein per gram liver (mg/g)") # File S2 MPPGL = 39.79067 (Simcyp v17 healthy volunteers)
    sf_liver <- fixed(65); label("Empirical liver clearance scaling factor, SF (unitless)") # Table S1 'SF 65, empirically determined'

    # Placental transfer (ex vivo cotyledon perfusion; Table 1)
    clapp_cot_mf <- fixed(1.03); label("Apparent ex vivo maternal-to-fetal cotyledon clearance (mL/min/cotyledon)") # Table 1 'CL app cot mf 1.03 +/- 0.06'
    clapp_cot_fm <- fixed(1.03); label("Apparent ex vivo fetal-to-maternal cotyledon clearance (mL/min/cotyledon)") # Table 1 'CL app cot fm 1.03 +/- 0.23'
    fu_perf_m <- fixed(0.048); label("Fraction unbound in maternal perfusion buffer (unitless)") # Table 1 'fu mf 0.048'
    fu_perf_f <- fixed(0.043); label("Fraction unbound in fetal perfusion buffer (unitless)") # Table 1 'fu fm 0.043'
    ncot <- fixed(30); label("Number of cotyledons per placenta (count)") # Table 1 'Ncot 30'

    # Amniotic-fluid exchange (Underwood 2005 flows scaled by A_SF)
    sf_amniotic <- fixed(3); label("Amniotic-fluid exchange scaling factor, A_SF (unitless)") # File S2 A_SF = 3 (fitted to one amniotic-fluid value, Discussion)
  })

  model({
    ka <- exp(lka)
    fdepot <- exp(lfdepot)
    kp_lung <- exp(lkp_lung)
    kp_adipose <- exp(lkp_adipose)
    kp_bone <- exp(lkp_bone)
    kp_brain <- exp(lkp_brain)
    kp_heart <- exp(lkp_heart)
    kp_kidney <- exp(lkp_kidney)
    kp_muscle <- exp(lkp_muscle)
    kp_skin <- exp(lkp_skin)
    kp_spleen <- exp(lkp_spleen)
    kp_gut <- exp(lkp_gut)
    kp_liver <- exp(lkp_liver)
    kp_rest <- exp(lkp_rest)
    kp_rest_fet <- exp(lkp_rest_fet)

    fub <- fu / bpr # File S2 fub_m = fup_m/BP
    fub_fet <- fu_fet / bpr_fet # File S2 fub_f = fup_f/BP_f

    # ---- Maternal physiology (Table S2; File S2 'POPULATION SPECIFIC') ----
    bw <- 61.1 + 0.24098 * GA + 0.0038 * GA^2 # body weight (kg); 73.7 at GA 34
    co <- 301 + 5.916 * GA - 0.088 * GA^2 # cardiac output (L/h); 400.4 at GA 34

    # Fractional tissue volumes (Gaohua 2012, File S2)
    v_lung <- 0.004 * bw
    v_adipose <- 0.385 * bw
    v_bone <- 0.025 * bw
    v_brain <- 0.016 * bw
    v_heart <- 0.003 * bw
    v_kidney <- 0.0037 * bw
    v_muscle <- 0.278 * bw
    v_skin <- 0.0298 * bw
    v_spleen <- 0.0018 * bw
    v_gut <- 0.015 * bw
    v_liver <- 0.02 * bw
    v_rest <- 0.1134 * bw
    v_arterial <- 0.0205 * bw
    v_venous <- 0.0411 * bw
    lw <- v_liver # liver weight (kg), assumed equal to liver volume

    # Fractional blood flows (Simcyp v17 healthy female, File S2)
    q_lung <- 1 * co
    q_adipose <- 0.085 * co
    q_bone <- 0.05 * co
    q_brain <- 0.12 * co
    q_heart <- 0.05 * co
    q_kidney <- 0.17 * co
    q_muscle <- 0.12 * co
    q_skin <- 0.05 * co
    q_spleen <- 0.03 * co
    q_ha <- 0.065 * co
    q_gut <- 0.17 * co
    q_liver <- q_spleen + q_gut + q_ha # total hepatic outflow
    q_rest <- 0.090 * co

    # ---- Fetal physiology (Table S2; File S2) ----
    v_tot_fet <- 0.01 * exp((0.955 / 0.0702) * (1 - exp(-0.0702 * GA))) / 1000 # fetal volume (L); 2.32 at GA 34
    v_blood_fet <- (11.2 * GA - 209.4) / 1000 # fetal blood (L); 0.171 at GA 34
    v_rest_fet <- v_tot_fet - v_blood_fet
    v_amniotic <- (1.9648 * GA - 1.2056 * GA^2 + 0.2064 * GA^3 - 0.0061 * GA^4 + 0.00005 * GA^5) / 1000 # amniotic fluid (L); 0.906 at GA 34
    co_fet <- 553 * v_tot_fet * (60 / 1000) # fetal cardiac output (L/h)
    fr_qpla_fet <- exp(3.35420863 + 0.000060601 * GA^3 - 0.000018693 * GA^3 * log(GA)) / 100 # fraction of fetal CO to placenta; 0.232 at GA 34
    q_rest_fet <- co_fet * (1 - fr_qpla_fet)

    # Amniotic-fluid exchange flows (Underwood 2005, L/day -> L/h)
    q_fet_amn <- (0.275 / sf_amniotic) / 24 # oral/nasal/tracheal/pulmonary secretion (kL)
    q_swallow <- (0.774 * sf_amniotic) / 24 # fetal swallowing (kSw)
    q_intramem <- (0.350 * sf_amniotic) / 24 # intramembranous uptake (kIntra)

    # ---- Hepatic clearance (File S2 'Hepatic clearance') ----
    enz_increase <- (100 + 2.9826 * GA - 0.0741 * GA^2) / 100 # CYP3A4 induction, UGT1A1 assumed equal
    clint_hep <- ((clint_ugt1a1 + clint_cyp3a4) / fumic) * mppgl * (lw * 1000) * (60 / 1e6) * enz_increase # L/h
    cl_liver <- clint_hep * sf_liver

    # ---- Placental clearance (Methods Eqs; File S2 'Cotyledon clearance') ----
    cl_plac_mf <- clapp_cot_mf / fu_perf_m * (60 / 1000) * fub * ncot # L/h
    cl_plac_fm <- clapp_cot_fm / fu_perf_f * (60 / 1000) * fub_fet * ncot # L/h

    # ---- Concentrations (mg/L) ----
    c_gut <- gut / v_gut
    c_lung <- lung / v_lung
    c_adipose <- adipose / v_adipose
    c_bone <- bone / v_bone
    c_brain <- brain / v_brain
    c_heart <- heart / v_heart
    c_kidney <- kidney / v_kidney
    c_muscle <- muscle / v_muscle
    c_skin <- skin / v_skin
    c_spleen <- spleen / v_spleen
    c_liver <- liver / v_liver
    c_rest <- rest / v_rest
    c_art <- arterial / v_arterial
    c_ven <- venous / v_venous
    c_blood_fet <- blood_fet / v_blood_fet
    c_rest_fet <- rest_fet / v_rest_fet
    c_amniotic <- amniotic / v_amniotic
    cu_liver <- c_liver * fu # File S2 Cu_liv = Cliver * fup_m

    # ---- Maternal ODEs (File S2 'DIFFERENTIAL EQUATIONS', mg/h) ----
    d/dt(depot) <- -ka * depot
    f(depot) <- fdepot
    d/dt(gut) <- ka * depot + q_gut * (c_art - c_gut / kp_gut * bpr)
    d/dt(lung) <- q_lung * c_ven - q_lung * (c_lung / kp_lung * bpr)
    d/dt(adipose) <- q_adipose * (c_art - c_adipose / kp_adipose * bpr)
    d/dt(bone) <- q_bone * (c_art - c_bone / kp_bone * bpr)
    d/dt(brain) <- q_brain * (c_art - c_brain / kp_brain * bpr)
    d/dt(heart) <- q_heart * (c_art - c_heart / kp_heart * bpr)
    d/dt(kidney) <- q_kidney * (c_art - c_kidney / kp_kidney * bpr)
    d/dt(muscle) <- q_muscle * (c_art - c_muscle / kp_muscle * bpr)
    d/dt(skin) <- q_skin * (c_art - c_skin / kp_skin * bpr)
    d/dt(spleen) <- q_spleen * (c_art - c_spleen / kp_spleen * bpr)
    d/dt(liver) <- q_ha * c_art + q_gut * (c_gut / kp_gut * bpr) + q_spleen * (c_spleen / kp_spleen * bpr) -
      q_liver * (c_liver / kp_liver * bpr) - cl_liver * cu_liver
    d/dt(rest) <- q_rest * (c_art - c_rest / kp_rest * bpr)
    d/dt(arterial) <- q_lung * (c_lung / kp_lung * bpr) - q_lung * c_art - cl_plac_mf * c_art
    d/dt(venous) <- q_adipose * (c_adipose / kp_adipose * bpr) + q_bone * (c_bone / kp_bone * bpr) +
      q_brain * (c_brain / kp_brain * bpr) + q_heart * (c_heart / kp_heart * bpr) +
      q_kidney * (c_kidney / kp_kidney * bpr) + q_muscle * (c_muscle / kp_muscle * bpr) +
      q_skin * (c_skin / kp_skin * bpr) + q_liver * (c_liver / kp_liver * bpr) +
      q_rest * (c_rest / kp_rest * bpr) + cl_plac_fm * c_blood_fet - q_lung * c_ven

    # ---- Fetoplacental ODEs (File S2) ----
    d/dt(blood_fet) <- cl_plac_mf * c_art + q_rest_fet * (c_rest_fet / kp_rest_fet * bpr_fet - c_blood_fet) +
      c_amniotic * q_intramem + c_amniotic * q_swallow - c_blood_fet * q_fet_amn - cl_plac_fm * c_blood_fet
    d/dt(rest_fet) <- q_rest_fet * (c_blood_fet - c_rest_fet / kp_rest_fet * bpr_fet)
    d/dt(amniotic) <- c_blood_fet * q_fet_amn - c_amniotic * q_intramem - c_amniotic * q_swallow

    # ---- Outputs ----
    Cc <- c_ven / bpr # maternal venous plasma (File S2 Cven_plasmaM)
    Cplasma_fet <- c_blood_fet / bpr_fet # fetal plasma (File S2 CplasmaF)
    Camniotic <- c_amniotic # amniotic fluid (File S2 Camf)
  })
}
