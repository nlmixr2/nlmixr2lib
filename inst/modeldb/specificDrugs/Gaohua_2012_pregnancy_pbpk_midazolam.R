Gaohua_2012_pregnancy_pbpk_midazolam <- function() {
  description <- paste(
    "PBPK (whole-body, 14-compartment pregnancy p-PBPK adapted from the",
    "Simcyp Simulator version 11 full-PBPK platform). Midazolam (CYP3A4",
    "substrate) disposition in healthy Caucasian women aged 20-40 years,",
    "with gestational-age-dependent maternal physiology and hepatic CYP3A4",
    "activity. The 14 ODE compartments are arterial blood, venous blood,",
    "lung, adipose, bone, brain, heart, kidney, gut, liver, muscle, skin,",
    "spleen, and a lumped fetoplacental unit (fetus + placenta + amniotic",
    "fluid + membranes + umbilical cord) per Gaohua 2012 Figure 1; the",
    "uterus and mammary glands are merged into the muscle compartment, so",
    "muscle volume and flow are computed as the residual that balances total",
    "body weight and cardiac output during pregnancy. Time-varying physiology",
    "(cardiac output, body weight, plasma / RBC volumes, hematocrit, serum",
    "albumin, skin / adipose / renal / fetoplacental blood flows, and CYP1A2",
    "/ CYP2D6 / CYP3A4 enzyme activities) follows the polynomial formula",
    "X = X0 * (a0 + a1*GA + a2*GA^2 + a3*GA^3) in Table 2; the fetoplacental",
    "volume follows the Gompertz curve in Eq. 1. Drug-specific values for",
    "fa, Fg, ka, fu, B:P, basal CL_int,H, the CYP fractional contributions",
    "A_1A2 / A_2D6 / A_3A4, and the 12 tissue:plasma partition coefficients",
    "(Rodgers and Rowland) are from Tables 3-4. Set covariate GA = 0 to",
    "simulate the non-pregnant reference woman; values 0 < GA < 40 simulate",
    "any gestational stage. The model is a perfusion-limited typical-value",
    "PBPK forward simulation; the paper added no IIV or residual-error model."
  )
  reference <- paste(
    "Gaohua L, Abduljalil K, Jamei M, Johnson TN, Rostami-Hodjegan A.",
    "A pregnancy physiologically based pharmacokinetic (p-PBPK) model for",
    "disposition of drugs metabolized by CYP1A2, CYP2D6 and CYP3A4.",
    "Br J Clin Pharmacol. 2012;74(5):873-885.",
    "doi:10.1111/j.1365-2125.2012.04363.x.",
    sep = " "
  )
  vignette <- "Gaohua_2012_pregnancy_pbpk"
  units    <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  paper_specific_compartments <- c("preg")

  covariateData <- list(
    GA = list(
      description        = "Gestational age at the time of dose administration",
      units              = "weeks",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Gestational age (in weeks) of the mother when the simulation starts.",
        "Inside model() the running variable gaweek evolves as",
        "gaweek = GA + t/(24*7) so physiological parameters drift over a long",
        "simulation, but for the typical single-dose PK time window (hours to",
        "a few days) gaweek is effectively constant at GA. Set GA = 0 to",
        "recover the non-pregnant",
        "reference physiology (Q_preg = 0, V_preg = 0.01 L); GA in the",
        "1-40 weeks range simulates the corresponding pregnancy stage.",
        "Per the paper Methods (Eq. 3): t in hours, 24*7 = 168 hours per week."
      ),
      source_name        = "GA"
    )
  )

  population <- list(
    species        = "human (pregnant women)",
    n_subjects     = NA_integer_,
    n_studies      = 3L,
    age_range      = "20-40 years (basal reference); validated against pregnancy stages 28-37 weeks",
    weight_range   = "60 kg basal non-pregnant; rises to ~76 kg at 36 weeks per Table 2",
    sex_female_pct = 100,
    race_ethnicity = c(White = 100),
    disease_state  = paste(
      "Healthy pregnant Caucasian women. The validation cohorts in the Gaohua",
      "2012 paper were caffeine (n = 8 pregnant at 36 +/- 3 weeks gestation,",
      "n = 4 nonpregnant; Brazier 1983 ref 24), metoprolol (n = 5 pregnant at",
      "~37 weeks, then 3-6 months postpartum; Hogstedt 1985 ref 28), and",
      "midazolam (n = 13 pregnant at 28-32 weeks then 6-10 weeks postpartum;",
      "Hebert 2008 ref 29). All three drugs were single oral or IV doses; this",
      "model file extracts the midazolam case and is paired with",
      "Gaohua_2012_pregnancy_pbpk_caffeine and",
      "Gaohua_2012_pregnancy_pbpk_metoprolol."
    ),
    dose_range     = "Midazolam 2 mg oral single dose (Hebert 2008 protocol; Gaohua 2012 Fig 4, Table 5)",
    regions        = "Europe (UK / France)",
    notes          = paste(
      "Basal physiology was taken from Simcyp Simulator version 11 for the",
      "healthy non-pregnant 20-40 year old Caucasian woman cohort. The 14-",
      "compartment perfusion-limited Simcyp full-PBPK structure (Jamei 2009",
      "ref 22) was extended with one extra perfusion-limited compartment",
      "representing the fetoplacental unit, time-varying physiology per the",
      "Abduljalil 2012 meta-analysis (ref 3), and gestational-age-dependent",
      "hepatic CYP1A2 / CYP2D6 / CYP3A4 activity. Implementation: Matlab",
      "Simulink R2010a in the source paper; this file ports the same 14",
      "ODEs and parameter tables to nlmixr2 / rxode2. No individual-level",
      "variability is reported; the paper states 'Unlike the Simcyp Simulator,",
      "no variability terms were added to the model.'"
    )
  )

  ini({
    # All structural parameters are FIXED at the Gaohua 2012 paper values
    # (Simcyp library entries for the three case drugs and the Rodgers and
    # Rowland tissue:plasma partition coefficients). The propSd residual
    # error is a placeholder for syntactic completeness; see Assumptions.

    # Absorption parameters (Gaohua 2012 Table 3, midazolam column).
    lka       <- fixed(log(3.04));     label("First-order absorption rate constant (1/h)")                 # Gaohua 2012 Table 3 (midazolam ka = 3.04 1/h)
    fa        <- fixed(0.88);          label("Fraction available for absorption from dosage form")         # Gaohua 2012 Table 3 (midazolam fa = 0.88)
    fg        <- fixed(0.59);          label("Gut availability (fraction surviving gut wall first pass)")  # Gaohua 2012 Table 3 (midazolam Fg = 0.59)

    # Plasma protein binding and blood:plasma ratio at basal Hct = 38.3% and
    # serum albumin = 44.9 g/L (Gaohua 2012 Tables 2-3). Held constant over
    # GA in this implementation; see Assumptions in the vignette.
    fu        <- fixed(0.032);         label("Fraction unbound in plasma (basal, held constant over GA)")  # Gaohua 2012 Table 3 (midazolam fu = 0.032)
    bp        <- fixed(0.664);         label("Blood:plasma concentration ratio (basal, held constant)")    # Gaohua 2012 Table 3 (midazolam B:P = 0.664)

    # Hepatic intrinsic clearance at basal CYP activity, log-transformed so
    # that any user re-fit retains positivity (paper value is fixed here).
    lcl_int_h <- fixed(log(1583));     label("Basal hepatic intrinsic clearance CL_int,H0 (L/h)")          # Gaohua 2012 Table 3 (midazolam CL_int,H = 1583 L/h)

    # Fractional contribution of each CYP isoform to hepatic intrinsic
    # clearance (Gaohua 2012 Eq. 6; A_1A2 + A_2D6 + A_3A4 = 1 by Table 3).
    a_cyp1a2  <- fixed(0.0);           label("Fractional contribution of CYP1A2 to CL_int,H")              # Gaohua 2012 Table 3 (midazolam A_1A2 = 0)
    a_cyp2d6  <- fixed(0.0);           label("Fractional contribution of CYP2D6 to CL_int,H")              # Gaohua 2012 Table 3 (midazolam A_2D6 = 0)
    a_cyp3a4  <- fixed(1.0);           label("Fractional contribution of CYP3A4 to CL_int,H")              # Gaohua 2012 Table 3 (midazolam A_3A4 = 1)

    # Tissue:plasma partition coefficients (Rodgers and Rowland, refs 31-32).
    # Held constant across GA; the fetoplacental Kp is assumed equal to brain
    # per Gaohua 2012 Methods (page 877, paragraph after Table 4).
    kp_adipose <- fixed(9.334);   label("Tissue:plasma partition coefficient, adipose")     # Gaohua 2012 Table 4 (midazolam adipose 9.334)
    kp_bone    <- fixed(7.712);   label("Tissue:plasma partition coefficient, bone")        # Gaohua 2012 Table 4 (midazolam bone 7.712)
    kp_brain   <- fixed(7.04);    label("Tissue:plasma partition coefficient, brain")       # Gaohua 2012 Table 4 (midazolam brain 7.04)
    kp_gut     <- fixed(5.62);    label("Tissue:plasma partition coefficient, gut")         # Gaohua 2012 Table 4 (midazolam gut 5.62)
    kp_heart   <- fixed(1.805);   label("Tissue:plasma partition coefficient, heart")       # Gaohua 2012 Table 4 (midazolam heart 1.805)
    kp_kidney  <- fixed(2.725);   label("Tissue:plasma partition coefficient, kidney")      # Gaohua 2012 Table 4 (midazolam kidney 2.725)
    kp_liver   <- fixed(3.4);     label("Tissue:plasma partition coefficient, liver")       # Gaohua 2012 Table 4 (midazolam liver 3.4)
    kp_lung    <- fixed(0.728);   label("Tissue:plasma partition coefficient, lung")        # Gaohua 2012 Table 4 (midazolam lung 0.728)
    kp_muscle  <- fixed(2.842);   label("Tissue:plasma partition coefficient, muscle")      # Gaohua 2012 Table 4 (midazolam muscle 2.842)
    kp_skin    <- fixed(3.436);   label("Tissue:plasma partition coefficient, skin")        # Gaohua 2012 Table 4 (midazolam skin 3.436)
    kp_spleen  <- fixed(2.757);   label("Tissue:plasma partition coefficient, spleen")      # Gaohua 2012 Table 4 (midazolam spleen 2.757)
    kp_preg    <- fixed(7.04);    label("Tissue:plasma partition coefficient, fetoplacental unit") # Gaohua 2012 Table 4 footnote (assumed = brain)

    # Placeholder residual error for the typical-value forward simulation.
    # The Gaohua 2012 paper did not fit an error model; this is provided so
    # rxode2 stochastic simulations are syntactically valid.
    propSd <- fixed(0.10); label("Proportional residual error placeholder (fraction)")  # not fit; placeholder for forward simulation
  })

  model({
    # ===== Gestational age (weeks) at current PK time t (hours) =====
    # Gaohua 2012 Eq. 3: gaweek(t) = GA + t / (24*7). The covariate GA is
    # the mother's gestational age at the time of dose administration; the
    # running value gaweek evolves during the simulation. The pregnancy time
    # scale (weeks) is long compared with the PK time scale (hours), so
    # gaweek is effectively constant at GA over a typical single-dose
    # simulation, but the formula keeps the model consistent for multi-day
    # or multi-week sims.
    gaweek <- GA + t / 168.0

    # ===== Time-varying maternal physiology (Gaohua 2012 Tables 1-2, Eq. 2)
    # Polynomial form X = X0 * (a0 + a1*gaweek + a2*gaweek^2 + a3*gaweek^3);
    # a0 = 1 for all but Q_preg (a0 = 0 so flow is zero at conception).

    Qc        <- 300   * (1 + 0.01965*gaweek - 0.000292*gaweek*gaweek)                                       # Table 2 row 1 (cardiac output, L/h)
    BW        <- 60    * (1 + 0.00560*gaweek + 0.000054*gaweek*gaweek)                                       # Table 2 row 2 (total bodyweight, kg)
    V_adipose <- 25.89 * (1 + 0.00035*gaweek + 0.000152*gaweek*gaweek)                                       # Table 2 row 3 (total fat mass V_adip, kg)
    V_plas    <- 2.71  * (1 - 0.00892*gaweek + 0.00168*gaweek*gaweek - 0.000028*gaweek*gaweek*gaweek)        # Table 2 row 4 (plasma volume, L)
    V_rbc     <- 1.45  * (1 + 0.00658*gaweek)                                                                # Table 2 row 5 (RBC volume, L)
    Q_skin    <- 14.975* (1 + 0.02882*gaweek)                                                                # Table 2 row 8 (skin blood flow, L/h)
    Q_adipose <- 25.458* (1 + 0.01542*gaweek - 0.000220*gaweek*gaweek)                                       # Table 2 row 9 (adipose blood flow, L/h)
    Q_kidney  <- 50.915* (1 + 0.05022*gaweek - 0.00125*gaweek*gaweek)                                        # Table 2 row 10 (renal blood flow, L/h)
    # Table 2 row 11 (fetoplacental blood flow, L/h; a0 = 0). The polynomial
    # is negative for gaweek in roughly (0, 4) weeks (Q_preg ~ -0.29 at
    # gaweek = 1), an artefact of the cubic fit calibrated against term-
    # gestation data. Clamp at zero so the non-pregnant / first-trimester
    # window has Q_preg = 0 (well-stirred fetoplacental ODE collapses to
    # no transfer).
    Q_preg_raw <- -0.4051*gaweek + 0.1188*gaweek*gaweek - 0.0019*gaweek*gaweek*gaweek
    Q_preg     <- max(0, Q_preg_raw)

    # Gompertz growth curve for the fetoplacental compartment volume
    # (Gaohua 2012 Eq. 1): V_preg = a * exp((b/c) * (1 - exp(-c*gaweek))).
    # Constants a = 0.01, b = 0.37, c = 0.052 per the paragraph below Eq. 1.
    V_preg <- 0.01 * exp((0.37 / 0.052) * (1 - exp(-0.052 * gaweek)))

    # CYP enzyme activity scaling factors (Table 2 rows 14-16). At gaweek = 0
    # each factor equals 1.0, recovering the basal CL_int,H.
    alpha_cyp1a2 <- 1 + (-0.03581)*gaweek + 0.00050 *gaweek*gaweek
    alpha_cyp2d6 <- 1 +   0.02270 *gaweek + (-0.00035)*gaweek*gaweek
    alpha_cyp3a4 <- 1 +   0.02983 *gaweek + (-0.00074)*gaweek*gaweek

    # ===== Compartment volumes (L) - Table 1 fractional values for basal
    # static tissues, scaled to the time-varying total body weight. =====
    V_bone   <- 0.025  * BW
    V_brain  <- 0.016  * BW
    V_heart  <- 0.003  * BW
    V_kidney <- 0.0037 * BW
    V_liver  <- 0.02   * BW
    V_lung   <- 0.004  * BW
    V_skin   <- 0.0298 * BW
    V_spleen <- 0.0018 * BW
    V_gut    <- 0.015  * BW

    # Arterial:venous = 1:2 split of total blood volume (paper page 875,
    # paragraph "the total blood volume ... was further divided into the
    # arterial and venous parts in the proportion of 1:2").
    V_blood    <- V_plas + V_rbc
    V_arterial <- V_blood / 3
    V_venous   <- V_blood * 2 / 3

    # Muscle volume is the balancing remainder (paper page 875: "the tissue
    # volume and blood flow rate of the muscle compartment were used to
    # balance the allocation of total bodyweight and cardiac output").
    V_muscle <- BW - (V_adipose + V_bone + V_brain + V_heart + V_kidney +
                      V_liver + V_lung + V_skin + V_spleen + V_gut +
                      V_arterial + V_venous + V_preg)

    # ===== Compartment blood flows (L/h) - Table 1 fractions scaled to the
    # time-varying cardiac output. Hepatic artery = Q_liver,total - Q_gut -
    # Q_spleen = (0.28 - 0.17 - 0.05) * Qc = 0.06 * Qc. =====
    Q_bone         <- 0.05 * Qc
    Q_brain        <- 0.12 * Qc
    Q_heart        <- 0.05 * Qc
    Q_gut          <- 0.17 * Qc
    Q_spleen       <- 0.05 * Qc
    Q_HA           <- 0.06 * Qc
    Q_liver_total  <- Q_HA + Q_gut + Q_spleen
    Q_lung         <- Qc

    # Muscle blood flow is the balancing remainder for cardiac output.
    Q_muscle <- Qc - (Q_adipose + Q_bone + Q_brain + Q_heart + Q_kidney +
                      Q_skin + Q_spleen + Q_gut + Q_HA + Q_preg)

    # ===== Individual structural parameters =====
    ka          <- exp(lka)
    cl_int_h0   <- exp(lcl_int_h)
    cl_int_h    <- cl_int_h0 * (alpha_cyp1a2 * a_cyp1a2 +
                                alpha_cyp2d6 * a_cyp2d6 +
                                alpha_cyp3a4 * a_cyp3a4)

    # ===== Tissue and blood concentrations =====
    # State variables hold AMOUNT of drug (dose units, mg) in each compartment.
    # C_t = mass / volume = mg/L. Blood concentration leaving tissue t is
    # C_blood_out_t = C_t * bp / Kp_t per Gaohua 2012 Eq. 5 rearrangement
    # (paper writes C_preg / (K_PregP / B:P) = C_preg * B:P / K_PregP).
    C_adipose <- adipose / V_adipose
    C_bone    <- bone    / V_bone
    C_brain   <- brain   / V_brain
    C_heart   <- heart   / V_heart
    C_kidney  <- kidney  / V_kidney
    C_gut     <- gut     / V_gut
    C_liver   <- liver   / V_liver
    C_lung    <- lung    / V_lung
    C_muscle  <- muscle  / V_muscle
    C_skin    <- skin    / V_skin
    C_spleen  <- spleen  / V_spleen
    C_preg    <- preg    / V_preg
    C_art     <- arterial / V_arterial
    C_ven     <- venous   / V_venous

    Cb_adipose <- C_adipose * bp / kp_adipose
    Cb_bone    <- C_bone    * bp / kp_bone
    Cb_brain   <- C_brain   * bp / kp_brain
    Cb_heart   <- C_heart   * bp / kp_heart
    Cb_kidney  <- C_kidney  * bp / kp_kidney
    Cb_gut     <- C_gut     * bp / kp_gut
    Cb_liver   <- C_liver   * bp / kp_liver
    Cb_lung    <- C_lung    * bp / kp_lung
    Cb_muscle  <- C_muscle  * bp / kp_muscle
    Cb_skin    <- C_skin    * bp / kp_skin
    Cb_spleen  <- C_spleen  * bp / kp_spleen
    Cb_preg    <- C_preg    * bp / kp_preg

    # ===== ODE system (Gaohua 2012 Eqs. 4-6 generalised across the 14
    # maternal + 1 fetoplacental compartments; depot adds oral absorption) =====

    # Oral absorption: depot empties at rate ka; fraction fa * fg reaches
    # the gut tissue (the residual fa*(1-fg) is lost to gut wall metabolism).
    d/dt(depot) <- -ka * depot

    # Lung sits between venous and arterial blood: receives full cardiac
    # output from venous, returns to arterial after one well-stirred pass.
    d/dt(lung)  <- Q_lung * (C_ven - Cb_lung)

    # Non-portal peripheral tissues drain to venous blood.
    d/dt(adipose) <- Q_adipose * (C_art - Cb_adipose)
    d/dt(bone)    <- Q_bone    * (C_art - Cb_bone)
    d/dt(brain)   <- Q_brain   * (C_art - Cb_brain)
    d/dt(heart)   <- Q_heart   * (C_art - Cb_heart)
    d/dt(kidney)  <- Q_kidney  * (C_art - Cb_kidney)

    # Gut: arterial input plus oral absorption (fa * fg surviving gut wall);
    # drains to liver via portal vein.
    d/dt(gut) <- ka * depot * fa * fg + Q_gut * (C_art - Cb_gut)

    # Liver: hepatic artery + portal (gut + spleen) inflow; hepatic vein
    # outflow back to venous; intrinsic-clearance loss on unbound plasma
    # concentration C_liver / kp_liver.
    d/dt(liver) <- Q_HA      * C_art      +
                   Q_gut     * Cb_gut     +
                   Q_spleen  * Cb_spleen  -
                   Q_liver_total * Cb_liver -
                   cl_int_h * fu * (C_liver / kp_liver)

    d/dt(muscle) <- Q_muscle * (C_art - Cb_muscle)
    d/dt(skin)   <- Q_skin   * (C_art - Cb_skin)

    # Spleen: drains to liver via portal vein.
    d/dt(spleen) <- Q_spleen * (C_art - Cb_spleen)

    # Fetoplacental unit (Gaohua 2012 Eq. 5): perfusion-limited compartment
    # in parallel with the other maternal tissues.
    d/dt(preg) <- Q_preg * (C_art - Cb_preg)

    # Arterial blood: receives oxygenated blood from lung at rate Q_lung;
    # drains to all peripheral tissues plus the hepatic artery. By mass
    # balance with Q_muscle as the residual, total arterial outflow = Qc.
    d/dt(arterial) <- Q_lung * (Cb_lung - C_art)

    # Venous blood: receives non-portal tissue outflows plus the combined
    # hepatic vein flow (Q_liver_total * Cb_liver) and returns to lung.
    d/dt(venous) <- Q_adipose * Cb_adipose +
                    Q_bone    * Cb_bone    +
                    Q_brain   * Cb_brain   +
                    Q_heart   * Cb_heart   +
                    Q_kidney  * Cb_kidney  +
                    Q_muscle  * Cb_muscle  +
                    Q_skin    * Cb_skin    +
                    Q_preg    * Cb_preg    +
                    Q_liver_total * Cb_liver -
                    Q_lung    * C_ven

    # Systemic plasma concentration observed in peripheral venous blood
    # (clinical sampling site). C_ven is the whole-blood concentration in
    # the venous pool; divide by B:P to get the plasma concentration the
    # paper reports against clinical observations.
    Cc <- C_ven / bp
    Cc ~ prop(propSd)
  })
}
