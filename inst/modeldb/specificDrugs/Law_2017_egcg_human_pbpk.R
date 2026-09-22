Law_2017_egcg_human_pbpk <- function() {
  description <- paste(
    "PBPK (whole-body, 13 blood-flow-limited tissue compartments).",
    "Single-catechin model for epigallocatechin gallate (EGCg) in humans (Law",
    "et al. 2017). Thirteen perfusion-limited compartments (lung, kidney,",
    "muscle, brain, liver, spleen, gut, bone, skin, heart, adipose, rest of",
    "body, blood split into a one-third arterial and a two-thirds venous",
    "pool) with first-order oral input into gut tissue after a lag time,",
    "renal elimination from the kidney, and an enterohepatic recycling loop",
    "in which a biliary clearance drains the liver through a",
    "three-sub-compartment bile duct into the gut lumen, from which drug is",
    "either reabsorbed into gut tissue or transported to faeces. Tissue",
    "volumes, blood flows and tissue/blood partition coefficients are the",
    "paper's reference human values; the model carries no IIV and is intended",
    "for typical-value simulation. Law et al. built their tea catechin",
    "MIXTURE model by linking three of these single-catechin models with no",
    "pharmacokinetic interaction between them, so the mixture is reproduced",
    "by simulating the sibling models together - see the vignette."
  )
  reference <- paste(
    "Law FCP, Yao M, Bi HC, Lam S. Physiologically based pharmacokinetic",
    "modeling of tea catechin mixture in rats and humans. Pharmacol Res",
    "Perspect. 2017;5(3):e00305. doi:10.1002/prp2.305."
  )
  vignette <- "Law_2017_teacatechins_pbpk"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # The bile duct is represented as a three-sub-compartment transit chain
  # (Appendix eq A5, n = 3, adapted from Bischoff et al. 1971 and Harrison and
  # Gibaldi 1977). Each state holds the amount in transit; dividing by the
  # residence time rt_bile gives the transfer rate R_j of the paper. No
  # canonical compartment family covers a biliary transit chain, so the states
  # are declared paper-specific.
  paper_specific_compartments <- c(
    "bile_transit1",
    "bile_transit2",
    "bile_transit3"
  )

  compartmentData <- list(
    depot = list(analyte = "EGCg", units = "mg", specimen = "administration site", verified = TRUE),
    arterial = list(analyte = "EGCg", units = "mg", specimen = "whole blood", verified = TRUE),
    venous = list(analyte = "EGCg", units = "mg", specimen = "whole blood", verified = TRUE),
    lung = list(analyte = "EGCg", units = "mg", specimen = "tissue", verified = TRUE),
    liver = list(analyte = "EGCg", units = "mg", specimen = "tissue", verified = TRUE),
    kidney = list(analyte = "EGCg", units = "mg", specimen = "tissue", verified = TRUE),
    gut = list(analyte = "EGCg", units = "mg", specimen = "tissue", verified = TRUE),
    gut_lumen = list(analyte = "EGCg", units = "mg", specimen = "faeces", verified = TRUE),
    bile_transit1 = list(analyte = "EGCg", units = "mg", specimen = "bile", verified = TRUE),
    bile_transit2 = list(analyte = "EGCg", units = "mg", specimen = "bile", verified = TRUE),
    bile_transit3 = list(analyte = "EGCg", units = "mg", specimen = "bile", verified = TRUE),
    brain = list(analyte = "EGCg", units = "mg", specimen = "tissue", verified = TRUE),
    heart = list(analyte = "EGCg", units = "mg", specimen = "tissue", verified = TRUE),
    spleen = list(analyte = "EGCg", units = "mg", specimen = "tissue", verified = TRUE),
    muscle = list(analyte = "EGCg", units = "mg", specimen = "tissue", verified = TRUE),
    skin = list(analyte = "EGCg", units = "mg", specimen = "tissue", verified = TRUE),
    bone = list(analyte = "EGCg", units = "mg", specimen = "tissue", verified = TRUE),
    adipose = list(analyte = "EGCg", units = "mg", specimen = "tissue", verified = TRUE),
    other = list(analyte = "EGCg", units = "mg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Body weight drives every physiological and pharmacokinetic scaling in",
        "the model. Tissue volumes are a fixed percentage of body weight and",
        "tissue blood flows a fixed percentage of cardiac output, which is itself",
        "allometric (Table 4 footnote 1 (CO = 16.1 * BW^0.75 L/h; 390 L/h at 70",
        "kg; Travis 1987)). Rate constants scale as k = kc * BW^-0.3 and",
        "clearances as CL = CLc * BW^0.66 (Table 5 footnotes 2-4; Travis 1987,",
        "Chiou et al. 1998). The gut lumen volume is a fixed absolute volume and",
        "does NOT scale with body weight. The simulated human studies used 72 kg",
        "(Chow et al. 2001, 2003; Lee et al. 2002 cohort 45-85 kg) and 75 kg",
        "(Chow et al. 2001 pure EGCg arm)."
      ),
      source_name = "BW"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 13L,
    n_studies = 4L,
    age_range = "adult volunteers (not further specified)",
    weight_range = "45-85 kg (simulations used 72 and 75 kg)",
    disease_state = paste(
      "Healthy adult volunteers. Law et al. did not generate new data; the",
      "human model was calibrated and validated against concentration-time",
      "profiles digitised from published studies - Chow et al. (2003), N = 8,",
      "72 kg, volunteers with Fitzpatrick type II or III skin; Chow et al.",
      "(2001), N = 5, 72-75 kg; Lee et al. (2002), 45-85 kg; and Chow et al.",
      "(2005). Cohort sizes for Lee et al. (2002) and Chow et al. (2005) are",
      "not reported in Law et al. (2017)."
    ),
    dose_range = paste(
      "Single oral doses. Pure EGCg 400 mg (Chow et al. 2001, 2003) or 2 mg/kg",
      "(Lee et al. 2002); Polyphenon E containing 400, 600, 800 or 1200 mg EGCg",
      "with EGC and EC in proportion; green tea solids 20 mg/kg containing EGCg",
      "13.9%, EGC 11.0% and EC 3.2% by weight (Lee et al. 2002)."
    ),
    regions = "not applicable (literature data reanalysis)",
    notes = paste(
      "Observed concentrations were read off the published figures with",
      "DigiMatic because the original data were unavailable. Goodness of fit",
      "was judged by mean absolute prediction error (MAPE) with a MAPE < 50%",
      "acceptance criterion; no residual-error model, no inter-individual",
      "variability and no parameter uncertainty are reported anywhere in the",
      "paper."
    )
  )

  ini({
    # Every value below is from Law et al. (2017) Table 5 (pharmacokinetic
    # parameters), Table 4 (tissue/blood partition coefficients) or Table 1
    # (blood/plasma ratio). Nothing was estimated here: Law et al. report no
    # standard errors, no confidence intervals and no random effects, so every
    # entry is fixed().

    # Table 5, square-bracket value (kac EGCg = 0.85 /h per kg^-0.3, Chow et al.
    # 2001); the brace value 1.1 simulates Chow et al. (2003) and the
    # parenthesis value 2.85 simulates Lee et al. (2002)
    lka <- fixed(log(0.85)); label("Absorption rate coefficient kac into gut tissue (/h per kg^-0.3)")
    # Table 5, square-bracket value (F EGCg = 0.12, Chow et al. 2001); 0.065
    # simulates Chow et al. (2003) and 0.07 simulates Lee et al. (2002)
    lfdepot <- fixed(log(0.12)); label("Empirical oral bioavailability factor F (unitless)")
    # Table 5 (tlag EGCg = 0.5 h)
    ltlag <- fixed(log(0.5)); label("Oral absorption lag time (h)")
    # Table 5 (Rt EGCg = 0.03 h)
    lrt_bile <- fixed(log(0.03)); label("Bile-duct sub-compartment residence time Rt (h)")
    # Table 5 (krac EGCg = 0.18)
    lkreab <- fixed(log(0.18)); label("Colonic reabsorption rate coefficient krac (/h per kg^-0.3)")
    # Table 5 (kfc EGCg = 25.6)
    lkfec <- fixed(log(25.6)); label("Faecal transport rate coefficient kfc (/h per kg^-0.3)")
    # Table 5 (CLbc EGCg = 2.7)
    lcl_nonren <- fixed(log(2.7)); label("Biliary clearance coefficient CLbc (L/h per kg^0.66)")
    # Table 5 (CLrc EGCg = 0.0023)
    lcl_renal <- fixed(log(0.0023)); label("Renal clearance coefficient CLrc (L/h per kg^0.66)")

    # Table 1 (BLPLR EGCg = 0.91; human values assumed equal to rat)
    bpr <- fixed(0.91); label("Blood/plasma concentration ratio BLPLR (unitless)")

    # Tissue/blood partition coefficients, predicted with the Poulin and Theil
    # tissue-composition model rather than measured (Table 4 footnote 3).
    # 'other' is the paper's 'rest of body' compartment.
    # Table 4 (adipose EGCg = 0.15)
    lkp_adipose <- fixed(log(0.15)); label("Tissue/blood partition coefficient, adipose (unitless)")
    # Table 4 (bone EGCg = 3.22)
    lkp_bone <- fixed(log(3.22)); label("Tissue/blood partition coefficient, bone (unitless)")
    # Table 4 (brain EGCg = 3.12)
    lkp_brain <- fixed(log(3.12)); label("Tissue/blood partition coefficient, brain (unitless)")
    # Table 4 (gut EGCg = 2.49)
    lkp_gut <- fixed(log(2.49)); label("Tissue/blood partition coefficient, gut (unitless)")
    # Table 4 (heart EGCg = 1)
    lkp_heart <- fixed(log(1)); label("Tissue/blood partition coefficient, heart (unitless)")
    # Table 4 (kidney EGCg = 1.38)
    lkp_kidney <- fixed(log(1.38)); label("Tissue/blood partition coefficient, kidney (unitless)")
    # Table 4 (liver EGCg = 2.05)
    lkp_liver <- fixed(log(2.05)); label("Tissue/blood partition coefficient, liver (unitless)")
    # Table 4 (lung EGCg = 0.57)
    lkp_lung <- fixed(log(0.57)); label("Tissue/blood partition coefficient, lung (unitless)")
    # Table 4 (muscle EGCg = 1.38)
    lkp_muscle <- fixed(log(1.38)); label("Tissue/blood partition coefficient, muscle (unitless)")
    # Table 4 (skin EGCg = 1.6)
    lkp_skin <- fixed(log(1.6)); label("Tissue/blood partition coefficient, skin (unitless)")
    # Table 4 (spleen EGCg = 1.4)
    lkp_spleen <- fixed(log(1.4)); label("Tissue/blood partition coefficient, spleen (unitless)")
    # Table 4 (other EGCg = 1)
    lkp_other <- fixed(log(1)); label("Tissue/blood partition coefficient, other (unitless)")

    # Law et al. report no residual-error model: goodness of fit was judged by
    # mean absolute prediction error (MAPE). nlmixr2 requires a residual term,
    # so propSd is fixed to the paper's own reported MAPE for this model. That
    # is a mean absolute relative deviation, not an estimated residual standard
    # deviation - do not read it as one. See the vignette Assumptions and
    # deviations section. Source: Results (13.2% MAPE, Figure 6, Chow et al.
    # 2003 calibration).
    propSd <- fixed(0.132); label("Proportional residual magnitude, from the reported MAPE (fraction)")
  })

  model({
    # Appendix of Law et al. (2017), equations A1-A10. States hold AMOUNTS in
    # mg; every concentration below is an amount divided by a volume in L, so
    # concentrations are mg/L (= ug/mL, the unit of the paper's figures).

    # Allometric exponents. Rate constants scale as k = kc * BW^-0.3 and
    # clearances as CL = CLc * BW^0.66 (Table 5 footnotes 2-4; Travis 1987,
    # Chiou et al. 1998). They are structural constants shared by every
    # parameter, so they live here rather than in ini().
    wt_exp_rate <- -0.3
    wt_exp_cl <- 0.66

    # Cardiac output, L/h. Table 4 footnote 1 (CO = 16.1 * BW^0.75 L/h; 390 L/h
    # at 70 kg; Travis 1987).
    co <- 16.1 * WT^0.75

    # Tissue volumes, L, as the paper's percentage of body weight (Table 4,
    # column 'Tissue volume (% BW)'). The 13 entries sum to 99.81% of body
    # weight.
    v_blood <- 0.0771 * WT
    v_adipose <- 0.12 * WT
    v_bone <- 0.0856 * WT
    v_brain <- 0.0002 * WT
    v_gut <- 0.0171 * WT
    v_heart <- 0.0047 * WT
    v_kidney <- 0.0044 * WT
    v_liver <- 0.0257 * WT
    v_lung <- 0.0076 * WT
    v_muscle <- 0.4 * WT
    v_skin <- 0.0371 * WT
    v_spleen <- 0.0026 * WT
    v_other <- 0.216 * WT

    # Total blood volume is split into a one-third arterial and a two-thirds
    # venous pool (Appendix, 'Blood compartment').
    v_arterial <- v_blood / 3
    v_venous <- 2 * v_blood / 3

    # Gut lumen volume, L. A fixed absolute volume that does not scale with body
    # weight (Table 4 footnote 2 (gut lumen volume 2.1 L; Bischoff et al.
    # 1971)).
    v_gut_lumen <- 2.1

    # Tissue blood flows, L/h, as the paper's percentage of cardiac output
    # (Table 4, column 'Blood flow (% CO)'). q_liver is TOTAL hepatic blood
    # flow; the hepatic arterial supply is q_liver - q_gut - q_spleen because
    # gut and spleen drain into the liver (Appendix eq A4). The nine flows that
    # return directly to venous blood plus q_liver sum to the cardiac output.
    q_adipose <- 0.05 * co
    q_bone <- 0.05 * co
    q_brain <- 0.12 * co
    q_gut <- 0.17 * co
    q_heart <- 0.04 * co
    q_kidney <- 0.19 * co
    q_liver <- 0.25 * co
    q_muscle <- 0.17 * co
    q_skin <- 0.05 * co
    q_spleen <- 0.02 * co
    q_other <- 0.08 * co

    # Individual parameters (no IIV: the paper is a typical-value model).
    ka <- exp(lka) * WT^wt_exp_rate
    kreab <- exp(lkreab) * WT^wt_exp_rate
    kfec <- exp(lkfec) * WT^wt_exp_rate
    cl_nonren <- exp(lcl_nonren) * WT^wt_exp_cl
    cl_renal <- exp(lcl_renal) * WT^wt_exp_cl
    fdepot <- exp(lfdepot)
    tlag <- exp(ltlag)
    rt_bile <- exp(lrt_bile)
    kp_adipose <- exp(lkp_adipose)
    kp_bone <- exp(lkp_bone)
    kp_brain <- exp(lkp_brain)
    kp_gut <- exp(lkp_gut)
    kp_heart <- exp(lkp_heart)
    kp_kidney <- exp(lkp_kidney)
    kp_liver <- exp(lkp_liver)
    kp_lung <- exp(lkp_lung)
    kp_muscle <- exp(lkp_muscle)
    kp_skin <- exp(lkp_skin)
    kp_spleen <- exp(lkp_spleen)
    kp_other <- exp(lkp_other)

    # Concentrations, mg/L.
    c_arterial <- arterial / v_arterial
    c_venous <- venous / v_venous
    c_adipose <- adipose / v_adipose
    c_bone <- bone / v_bone
    c_brain <- brain / v_brain
    c_gut <- gut / v_gut
    c_heart <- heart / v_heart
    c_kidney <- kidney / v_kidney
    c_liver <- liver / v_liver
    c_lung <- lung / v_lung
    c_muscle <- muscle / v_muscle
    c_skin <- skin / v_skin
    c_spleen <- spleen / v_spleen
    c_other <- other / v_other
    c_lung <- lung / v_lung
    c_gut_lumen <- gut_lumen / v_gut_lumen

    # Appendix eq A4: RAM = CLb * (C_liver / R_liver), the combined hepatic
    # metabolic and biliary secretion rate leaving the liver for the bile duct,
    # mg/h.
    ram <- cl_nonren * (c_liver / kp_liver)

    # Oral input. Appendix eq A7: RAO = ka * F * dose * exp(-ka * (t - tlag)),
    # an analytic first-order input into GUT TISSUE (not into the gut lumen). An
    # equivalent depot state is used here so rxode2 handles the dose record, the
    # bioavailability factor and the lag time; ka * depot is identical to the
    # paper's RAO.
    d/dt(depot) <- -ka * depot
    f(depot) <- fdepot
    alag(depot) <- tlag

    # Non-eliminating, blood-flow-limited organs. Appendix eq A1: V_X * dC_X/dt
    # = Q_X * (C_arterial - C_X / R_X).
    d/dt(adipose) <- q_adipose * (c_arterial - c_adipose / kp_adipose)
    d/dt(bone) <- q_bone * (c_arterial - c_bone / kp_bone)
    d/dt(brain) <- q_brain * (c_arterial - c_brain / kp_brain)
    d/dt(heart) <- q_heart * (c_arterial - c_heart / kp_heart)
    d/dt(muscle) <- q_muscle * (c_arterial - c_muscle / kp_muscle)
    d/dt(skin) <- q_skin * (c_arterial - c_skin / kp_skin)
    d/dt(spleen) <- q_spleen * (c_arterial - c_spleen / kp_spleen)
    d/dt(other) <- q_other * (c_arterial - c_other / kp_other)

    # Lung. Appendix eq A2: the lung sees the whole cardiac output arriving from
    # venous blood, and drains into arterial blood.
    d/dt(lung) <- co * (c_venous - c_lung / kp_lung)

    # Kidney, an eliminating organ. Appendix eq A3 adds the renal clearance
    # acting on the organ's blood-equivalent concentration C_kidney / R_kidney.
    d/dt(kidney) <- q_kidney * (c_arterial - c_kidney / kp_kidney) -
      cl_renal * (c_kidney / kp_kidney)

    # Liver, an eliminating organ. Appendix eq A4. Inflow is the hepatic artery
    # (q_liver - q_gut - q_spleen) plus the portal return from gut and spleen;
    # outflow is the whole hepatic flow plus the biliary flux RAM.
    d/dt(liver) <- (q_liver - q_gut - q_spleen) * c_arterial +
      q_gut * (c_gut / kp_gut) + q_spleen * (c_spleen / kp_spleen) -
      q_liver * (c_liver / kp_liver) - ram

    # Bile duct. Appendix eq A5: Rt * dR_j/dt = R_(j-1) - R_j with R_0 = RAM and
    # n = 3 sub-compartments. Written here on the AMOUNT scale, bile_transit_j =
    # Rt * R_j, which turns eq A5 into an ordinary transit chain with mean
    # residence time rt_bile per sub-compartment. The flux delivered to the gut
    # lumen is R_3 = bile_transit3 / rt_bile.
    d/dt(bile_transit1) <- ram - bile_transit1 / rt_bile
    d/dt(bile_transit2) <- (bile_transit1 - bile_transit2) / rt_bile
    d/dt(bile_transit3) <- (bile_transit2 - bile_transit3) / rt_bile

    # Gut lumen. Appendix eq A6: V_GC * dC_GC/dt = R_3 - kfec * C_GC * V_GT -
    # kreab * V_GC * C_GC. Note that the faecal-transport term is printed with
    # the GUT TISSUE volume V_GT, not the gut lumen volume V_GC; it is
    # transcribed here exactly as printed. Substituting V_GC changes the rat
    # terminal half-life by about a factor of two and the human profile by less
    # than 0.5%; see the vignette Errata.
    d/dt(gut_lumen) <- bile_transit3 / rt_bile -
      kfec * c_gut_lumen * v_gut - kreab * gut_lumen

    # Gut tissue. Appendix eq A7: perfusion term, plus reabsorption from the gut
    # lumen, plus the oral input RAO.
    d/dt(gut) <- q_gut * (c_arterial - c_gut / kp_gut) +
      kreab * gut_lumen + ka * depot

    # Arterial blood. Appendix eq A8: fed by the lung, drained to every tissue.
    d/dt(arterial) <- co * (c_lung / kp_lung - c_arterial)

    # Venous blood. Appendix eq A9: the venous return of every organ that drains
    # directly to the vena cava, plus the whole hepatic outflow (which already
    # carries the gut and spleen returns), minus the cardiac output leaving for
    # the lung.
    d/dt(venous) <- q_adipose * (c_adipose / kp_adipose) +
      q_bone * (c_bone / kp_bone) +
      q_brain * (c_brain / kp_brain) +
      q_heart * (c_heart / kp_heart) +
      q_kidney * (c_kidney / kp_kidney) +
      q_muscle * (c_muscle / kp_muscle) +
      q_skin * (c_skin / kp_skin) +
      q_liver * (c_liver / kp_liver) +
      q_other * (c_other / kp_other) -
      co * c_venous

    # Observation. Appendix eq A10: mixed venous PLASMA concentration is the
    # venous BLOOD concentration divided by the blood/plasma ratio. This is the
    # free (unconjugated) catechin concentration the paper's figures plot.
    Cc <- c_venous / bpr
    Cc ~ prop(propSd)
  })
}
