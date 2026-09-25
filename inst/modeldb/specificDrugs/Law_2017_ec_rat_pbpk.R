Law_2017_ec_rat_pbpk <- function() {
  description <- paste(
    "PBPK (whole-body, 13 blood-flow-limited tissue compartments).",
    "Preclinical (rat). Single-catechin model for epicatechin (EC) in rats",
    "(Law et al. 2017). Thirteen perfusion-limited compartments (lung,",
    "kidney, muscle, brain, liver, spleen, gut, bone, skin, heart, adipose,",
    "rest of body, blood split into a one-third arterial and a two-thirds",
    "venous pool) with first-order oral input into gut tissue after a lag",
    "time, renal elimination from the kidney, and an enterohepatic recycling",
    "loop in which a biliary clearance drains the liver through a",
    "three-sub-compartment bile duct into the gut lumen, from which drug is",
    "either reabsorbed into gut tissue or transported to faeces. Tissue",
    "volumes, blood flows and tissue/blood partition coefficients are the",
    "paper's reference rat values; the model carries no IIV and is intended",
    "for typical-value simulation. Law et al. built their tea catechin",
    "MIXTURE model by linking three of these single-catechin models with no",
    "pharmacokinetic interaction between them, so the mixture is reproduced",
    "by simulating the sibling models together - see the vignette. IMPORTANT:",
    "with Table 2 and Table 3 values exactly as printed this rat model does",
    "not quantitatively reproduce the paper's own Figures 4, 5 and 8; the",
    "human siblings do reproduce the paper's Table 7 predicted Cmax. See the",
    "vignette Errata."
  )
  reference <- paste(
    "Law FCP, Yao M, Bi HC, Lam S. Physiologically based pharmacokinetic",
    "modeling of tea catechin mixture in rats and humans. Pharmacol Res",
    "Perspect. 2017;5(3):e00305. doi:10.1002/prp2.305."
  )
  vignette <- "Law_2017_teacatechins_pbpk"
  units <- list(time = "min", dosing = "mg", concentration = "mg/L")

  # The bile duct is represented as a three-sub-compartment transit chain
  # (Appendix eq A5, n = 3, adapted from Bischoff et al. 1971 and Harrison and
  # Gibaldi 1977). Each state holds the amount in transit; dividing by the
  # per-sub-compartment residence time mtt_bile gives the transfer rate R_j of
  # the paper, so the total bile-duct delay is 3 * mtt_bile. bile_transit<n> is
  # a canonical chain family (operator ruling 2026-09-21, sidecar
  # oasweep_PMC5464336 q1); see inst/references/compartment-names.md.

  compartmentData <- list(
    depot = list(analyte = "EC", units = "mg", specimen = "administration site", verified = TRUE),
    arterial = list(analyte = "EC", units = "mg", specimen = "whole blood", verified = TRUE),
    venous = list(analyte = "EC", units = "mg", specimen = "whole blood", verified = TRUE),
    lung = list(analyte = "EC", units = "mg", specimen = "tissue", verified = TRUE),
    liver = list(analyte = "EC", units = "mg", specimen = "tissue", verified = TRUE),
    kidney = list(analyte = "EC", units = "mg", specimen = "tissue", verified = TRUE),
    gut = list(analyte = "EC", units = "mg", specimen = "tissue", verified = TRUE),
    gut_lumen = list(analyte = "EC", units = "mg", specimen = "faeces", verified = TRUE),
    bile_transit1 = list(analyte = "EC", units = "mg", specimen = "bile", verified = TRUE),
    bile_transit2 = list(analyte = "EC", units = "mg", specimen = "bile", verified = TRUE),
    bile_transit3 = list(analyte = "EC", units = "mg", specimen = "bile", verified = TRUE),
    brain = list(analyte = "EC", units = "mg", specimen = "tissue", verified = TRUE),
    heart = list(analyte = "EC", units = "mg", specimen = "tissue", verified = TRUE),
    spleen = list(analyte = "EC", units = "mg", specimen = "tissue", verified = TRUE),
    muscle = list(analyte = "EC", units = "mg", specimen = "tissue", verified = TRUE),
    skin = list(analyte = "EC", units = "mg", specimen = "tissue", verified = TRUE),
    bone = list(analyte = "EC", units = "mg", specimen = "tissue", verified = TRUE),
    adipose = list(analyte = "EC", units = "mg", specimen = "tissue", verified = TRUE),
    other = list(analyte = "EC", units = "mg", specimen = "tissue", verified = TRUE)
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
        "allometric (Table 2 footnote 1 (CO = 14.0 * BW^0.75 L/h; Travis 1987)).",
        "Rate constants scale as k = kc * BW^-0.3 and clearances as CL = CLc *",
        "BW^0.66 (Table 3 footnotes 2-4; Travis 1987, Chiou et al. 1998). The gut",
        "lumen volume is a fixed absolute volume and does NOT scale with body",
        "weight. Law et al. state a reference rat of 0.26 kg and scale to the",
        "mean weight of the rats in each simulated study (0.21-0.23 kg in Zhu et",
        "al. 2000, 0.31 kg in Chen et al. 1997)."
      ),
      source_name = "BW"
    )
  )

  population <- list(
    species = "rat (male Sprague-Dawley)",
    n_subjects = 6L,
    n_studies = 2L,
    age_range = "adult",
    weight_range = "0.21-0.31 kg",
    sex_female_pct = 0,
    disease_state = paste(
      "Healthy male Sprague-Dawley rats. Law et al. did not generate new data;",
      "the model was calibrated and validated against concentration-time",
      "profiles digitised from two published studies - Zhu et al. (2000), N = 6",
      "rats of 210-230 g with implanted jugular vein cannulae, and Chen et al.",
      "(1997), rats of 310 g sampled from the orbital sinus."
    ),
    dose_range = paste(
      "Single oral doses. Zhu et al. (2000): Polyphenon E containing EGCg 2500",
      "mg/kg, ECg 650 mg/kg and EC 250 mg/kg. Chen et al. (1997): pure EGCg 75",
      "mg/kg, or Polyphenon E containing EGCg 14.6 mg/kg, EGC 13.6 mg/kg and EC",
      "5.4 mg/kg."
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
    # Every value below is from Law et al. (2017) Table 3 (pharmacokinetic
    # parameters), Table 2 (tissue/blood partition coefficients) or Table 1
    # (blood/plasma ratio). Nothing was estimated here: Law et al. report no
    # standard errors, no confidence intervals and no random effects, so every
    # entry is fixed().
    #
    # IMPORTANT -- READ BEFORE USING THESE VALUES. With Table 2 and Table 3
    # exactly as printed, this RAT model does not quantitatively reproduce the
    # simulated curves Law et al. plot in their own Figures 4, 5 and 8 (EGCg
    # peaks ~2x low, ~3x late, and washes out ~4x slowly). The deviation is in
    # the SOURCE, not in the transcription or the implementation: every rat
    # number was re-read from a 300-dpi render of Tables 1-3, and the identical
    # ODE code reproduces the paper's human Table 7 predicted Cmax in 14 of 15
    # cells. Nothing here has been tuned to close the gap -- no single reading
    # of the printed parameters reconciles Cmax, Tmax and the terminal slope
    # together. Extracted with loud errata by operator ruling 2026-09-21
    # (sidecar oasweep_PMC5464336 q3 = A). See the vignette section 'The rat
    # models do not reproduce the paper's own rat figures'.

    # Table 3 (kac EC = 0.002 /min per kg^-0.3)
    lka <- fixed(log(0.002)); label("Absorption rate coefficient kac into gut tissue (/min per kg^-0.3)")
    # Table 3 (F EC = 0.13)
    lfdepot <- fixed(log(0.13)); label("Empirical oral bioavailability factor F (unitless)")
    # Table 3 (tlag EC = 5 min)
    ltlag <- fixed(log(5)); label("Oral absorption lag time (min)")
    # Table 3 (Rt EC = 2 min)
    lmtt_bile <- fixed(log(2)); label("Bile-duct sub-compartment residence time Rt (min)")
    # Table 3 (krac EC = 13.4)
    lkreab <- fixed(log(13.4)); label("Colonic reabsorption rate coefficient krac (/min per kg^-0.3)")
    # Table 3 (kfc EC = 0.13)
    lkfec <- fixed(log(0.13)); label("Faecal transport rate coefficient kfc (/min per kg^-0.3)")
    # Table 3 (CLbc EC = 8.7) (8.7 mL/min per kg^0.66 converted to L/min)
    lcl_nonren <- fixed(log(8.7 / 1000)); label("Biliary clearance coefficient CLbc (L/min per kg^0.66)")
    # Table 3 (CLrc EC = 4.5) (4.5 mL/min per kg^0.66 converted to L/min)
    lcl_renal <- fixed(log(4.5 / 1000)); label("Renal clearance coefficient CLrc (L/min per kg^0.66)")

    # Table 1 (BLPLR EC = 0.88, rat)
    bpr <- fixed(0.88); label("Blood/plasma concentration ratio BLPLR (unitless)")

    # Tissue/blood partition coefficients, predicted with the Poulin and Theil
    # tissue-composition model rather than measured (Table 2 footnote 3).
    # 'other' is the paper's 'rest of body' compartment.
    # Table 2 (adipose EC = 0.08)
    lkp_adipose <- fixed(log(0.08)); label("Tissue/blood partition coefficient, adipose (unitless)")
    # Table 2 (bone EC = 0.39)
    lkp_bone <- fixed(log(0.39)); label("Tissue/blood partition coefficient, bone (unitless)")
    # Table 2 (brain EC = 0.73)
    lkp_brain <- fixed(log(0.73)); label("Tissue/blood partition coefficient, brain (unitless)")
    # Table 2 (gut EC = 0.63)
    lkp_gut <- fixed(log(0.63)); label("Tissue/blood partition coefficient, gut (unitless)")
    # Table 2 (heart EC = 0.62)
    lkp_heart <- fixed(log(0.62)); label("Tissue/blood partition coefficient, heart (unitless)")
    # Table 2 (kidney EC = 0.63)
    lkp_kidney <- fixed(log(0.63)); label("Tissue/blood partition coefficient, kidney (unitless)")
    # Table 2 (liver EC = 0.59)
    lkp_liver <- fixed(log(0.59)); label("Tissue/blood partition coefficient, liver (unitless)")
    # Table 2 (lung EC = 0.65)
    lkp_lung <- fixed(log(0.65)); label("Tissue/blood partition coefficient, lung (unitless)")
    # Table 2 (muscle EC = 0.59)
    lkp_muscle <- fixed(log(0.59)); label("Tissue/blood partition coefficient, muscle (unitless)")
    # Table 2 (skin EC = 0.55)
    lkp_skin <- fixed(log(0.55)); label("Tissue/blood partition coefficient, skin (unitless)")
    # Table 2 (spleen EC = 0.63)
    lkp_spleen <- fixed(log(0.63)); label("Tissue/blood partition coefficient, spleen (unitless)")
    # Table 2 (other EC = 1)
    lkp_other <- fixed(log(1)); label("Tissue/blood partition coefficient, other (unitless)")

    # Law et al. report no residual-error model: goodness of fit was judged by
    # mean absolute prediction error (MAPE). nlmixr2 requires a residual term,
    # so propSd is fixed to the paper's own reported MAPE for this model. That
    # is a mean absolute relative deviation, not an estimated residual standard
    # deviation - do not read it as one. See the vignette Assumptions and
    # deviations section. Source: Results (33.9% overall MAPE, Figure 8, rat TCM
    # mixture).
    propSd <- fixed(0.339); label("Proportional residual magnitude, from the reported MAPE (fraction)")
  })

  model({
    # Appendix of Law et al. (2017), equations A1-A10. States hold AMOUNTS in
    # mg; every concentration below is an amount divided by a volume in L, so
    # concentrations are mg/L (= ug/mL, the unit of the paper's figures).

    # Allometric exponents. Rate constants scale as k = kc * BW^-0.3 and
    # clearances as CL = CLc * BW^0.66 (Table 3 footnotes 2-4; Travis 1987,
    # Chiou et al. 1998). They are structural constants shared by every
    # parameter, so they live here rather than in ini().
    wt_exp_rate <- -0.3
    wt_exp_cl <- 0.66

    # Cardiac output, L/min. Table 2 footnote 1 (CO = 14.0 * BW^0.75 L/h; Travis
    # 1987).
    co <- 14 * WT^0.75 / 60

    # Tissue volumes, L, as the paper's percentage of body weight (Table 2,
    # column 'Tissue volume (% BW)'). The 13 entries sum to 100% of body weight.
    v_blood <- 0.0816 * WT
    v_adipose <- 0.076 * WT
    v_bone <- 0.0415 * WT
    v_brain <- 0.0057 * WT
    v_gut <- 0.027 * WT
    v_heart <- 0.0033 * WT
    v_kidney <- 0.0073 * WT
    v_liver <- 0.0366 * WT
    v_lung <- 0.005 * WT
    v_muscle <- 0.404 * WT
    v_skin <- 0.19 * WT
    v_spleen <- 0.002 * WT
    v_other <- 0.12 * WT

    # Total blood volume is split into a one-third arterial and a two-thirds
    # venous pool (Appendix, 'Blood compartment').
    v_arterial <- v_blood / 3
    v_venous <- 2 * v_blood / 3

    # Gut lumen volume, L. A fixed absolute volume that does not scale with body
    # weight (Table 2 footnote 2 (gut lumen volume 0.0176 L; Angelo and
    # Pritchard 1987)).
    v_gut_lumen <- 0.0176

    # Tissue blood flows, L/min, as the paper's percentage of cardiac output
    # (Table 2, column 'Blood flow (% CO)'). q_liver is TOTAL hepatic blood
    # flow; the hepatic arterial supply is q_liver - q_gut - q_spleen because
    # gut and spleen drain into the liver (Appendix eq A4). The nine flows that
    # return directly to venous blood plus q_liver sum to the cardiac output.
    q_adipose <- 0.07 * co
    q_bone <- 0.122 * co
    q_brain <- 0.02 * co
    q_gut <- 0.131 * co
    q_heart <- 0.049 * co
    q_kidney <- 0.141 * co
    q_liver <- 0.175 * co
    q_muscle <- 0.278 * co
    q_skin <- 0.058 * co
    q_spleen <- 0.02 * co
    q_other <- 0.087 * co

    # Individual parameters (no IIV: the paper is a typical-value model).
    ka <- exp(lka) * WT^wt_exp_rate
    kreab <- exp(lkreab) * WT^wt_exp_rate
    kfec <- exp(lkfec) * WT^wt_exp_rate
    cl_nonren <- exp(lcl_nonren) * WT^wt_exp_cl
    cl_renal <- exp(lcl_renal) * WT^wt_exp_cl
    fdepot <- exp(lfdepot)
    tlag <- exp(ltlag)
    mtt_bile <- exp(lmtt_bile)
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
    # mg/min.
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
    # residence time mtt_bile per sub-compartment. The flux delivered to the gut
    # lumen is R_3 = bile_transit3 / mtt_bile.
    d/dt(bile_transit1) <- ram - bile_transit1 / mtt_bile
    d/dt(bile_transit2) <- (bile_transit1 - bile_transit2) / mtt_bile
    d/dt(bile_transit3) <- (bile_transit2 - bile_transit3) / mtt_bile

    # Gut lumen. Appendix eq A6: V_GC * dC_GC/dt = R_3 - kfec * C_GC * V_GT -
    # kreab * V_GC * C_GC. Note that the faecal-transport term is printed with
    # the GUT TISSUE volume V_GT, not the gut lumen volume V_GC; it is
    # transcribed here exactly as printed. Substituting V_GC changes the rat
    # terminal half-life by about a factor of two and the human profile by less
    # than 0.5%; see the vignette Errata.
    d/dt(gut_lumen) <- bile_transit3 / mtt_bile -
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
