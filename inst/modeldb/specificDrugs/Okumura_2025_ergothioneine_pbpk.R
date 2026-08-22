Okumura_2025_ergothioneine_pbpk <- function() {
  description <- paste(
    "PBPK (whole-body, OCTN1 transporter-mediated). Ergothioneine (EGT)",
    "disposition in healthy Japanese adults during repeated once-daily oral",
    "supplementation, used by Okumura 2025 to pick the lowest daily dose that",
    "reaches the plasma concentration associated with improved sleep quality.",
    "EGT is a dietary sulfur amino acid that mammals cannot synthesise; it",
    "crosses membranes only via its specific transporter OCTN1/SLC22A4, so",
    "every distribution and reabsorption step in the model is a saturable",
    "Michaelis-Menten process sharing one transporter Km. Plasma exchanges",
    "with three permeability-limited peripheral tissues (muscle, skin,",
    "adipose), each split into an extracellular and an intracellular space,",
    "and with a permeability-limited liver represented as five tandem",
    "extracellular/intracellular pairs to reproduce the near-flow-limited",
    "hepatic uptake seen with [3H]EGT in mice. Absorbed drug enters the first",
    "hepatic extracellular sub-compartment (portal delivery); gastrointestinal",
    "absorption is assumed complete. Elimination is renal: EGT is filtered at",
    "the glomerulus into a tubular duct, saturably reabsorbed back to plasma,",
    "and lost in urine. A red-blood-cell precursor pool takes EGT up from",
    "plasma saturably and matures into circulating red cells, which explains",
    "the long lag before whole-blood EGT rises. A continuous dietary intake",
    "rate feeds the gastrointestinal compartment alongside the supplement",
    "dose. The model is a fixed-effects least-squares fit to mean plasma and",
    "red-cell profiles (Napp software) with no between-subject variability;",
    "the values below are the final 'blush-up' fit of Supplementary Table S5,",
    "which added the 8 mg/day arm to the 5, 10 and 20 mg/day arms used for",
    "the Table 1 fit. Doses are supplied in umol (EGT MW 229.30 g/mol, so",
    "8 mg = 8 * 1000 / 229.30 = 34.89 umol). See the vignette Errata for the",
    "published discrepancies this extraction had to settle.",
    sep = " "
  )
  reference <- paste(
    "Okumura H, Araragi Y, Nishioka K, Yamashita R, Suzuki T, Watanabe H,",
    "Kato Y, Murayama N. Estimation and Validation of an Effective",
    "Ergothioneine Dose for Improved Sleep Quality Using Physiologically",
    "Based Pharmacokinetic Model. Food Sci Nutr. 2025;13(6):e70382.",
    "doi:10.1002/fsn3.70382.",
    "The ODE system is transcribed from the 'Actual model code used in Napp",
    "software' block of Supporting Information S1 (FSN3-13-e70382-s001.docx),",
    "cross-checked against the printed mass-balance equations in the same",
    "file; the initial-condition formulae are that file's 'Estimation of",
    "initial values' section; the fixed physiological constants are Table 1;",
    "and the estimated parameters are Supplementary Table S5.",
    "The five-tandem hepatic chain follows Watanabe T, Kusuhara H, Maeda K,",
    "Shitara Y, Sugiyama Y. Physiologically based pharmacokinetic modeling to",
    "predict transporter-mediated clearance and distribution of pravastatin in",
    "humans. J Pharmacol Exp Ther. 2009;328(2):652-62.",
    "doi:10.1124/jpet.108.146647, cited by Okumura 2025 as the source of that",
    "structure; the same chain is already used in this library by",
    "Aoki_2024_bosentan_pbpk.R.",
    sep = " "
  )
  vignette <- "Okumura_2025_ergothioneine_pbpk"

  # The five tandem hepatic extracellular / intracellular pairs are the
  # Watanabe 2009 dispersion idiom rather than a canonical single-organ
  # state, and are declared here exactly as the sibling model
  # Aoki_2024_bosentan_pbpk.R declares them. `rbc_precursor` is the
  # erythroid-precursor EGT pool in the bone marrow (an AMOUNT, umol) and
  # `renal_duct` is the proximal-tubular lumen; neither has a canonical
  # entry yet. Note the deliberate distinction between `rbc_precursor`
  # (amount of EGT in precursor cells) and `rbc_egt` (concentration of EGT
  # inside mature circulating red cells), which IS a canonical member of
  # the `rbc_<analyte>` family.
  paper_specific_compartments <- c(
    "is_liver1", "is_liver2", "is_liver3", "is_liver4", "is_liver5",
    "int_liver1", "int_liver2", "int_liver3", "int_liver4", "int_liver5",
    "rbc_precursor", "renal_duct"
  )

  # Time in hours (the Napp code steps in hours and doses every 24). The
  # model works in molar units throughout: tissue states hold
  # concentrations in umol/L (micromolar, the unit every EGT concentration
  # is reported in) and the two amount states hold umol, so doses must be
  # supplied in umol.
  units <- list(time = "h", dosing = "umol", concentration = "umol/L")

  # No covariates. Every physiological volume and flow in Table 1 is a
  # fixed constant for one typical 70 kg adult (Supplementary Table S5
  # footnote a: "Each value represents per human value assuming 70 kg body
  # weight"); the model carries no body-weight, age or sex term.
  covariateData <- list()

  # Screened in the paper's baseline tables but never entered into the
  # model, so documented here rather than in covariateData.
  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = paste(
        "Reported as a baseline characteristic (Study 1 Supplementary",
        "Table S3: 61.2 +/- 10.7 and 65.7 +/- 12.5 kg; Study 2 Table 3:",
        "61.3 +/- 11.1 and 61.6 +/- 12.5 kg) and used only to state that",
        "the physiological constants are for a 70 kg adult. No weight",
        "term was estimated."
      )
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = paste(
        "Baseline characteristic only (Study 1: 42.4 +/- 7.9 and 42.1",
        "+/- 8.1 years; Study 2: 53.5 +/- 6.3 and 53.4 +/- 6.0 years).",
        "The Discussion notes that plasma EGT is lower in older people",
        "(Cheah 2016) but no age term was estimated."
      )
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "unitless",
      type = "categorical",
      notes = paste(
        "Balanced by design at randomisation (Study 1: 8/28 female;",
        "Study 2: 51/100 female) but not entered into the model."
      )
    )
  )

  compartmentData <- list(
    depot         = list(analyte = "ergothioneine", units = "umol", specimen = "administration site", verified = TRUE),
    central       = list(analyte = "ergothioneine", units = "umol/L", specimen = "plasma", verified = TRUE),
    is_liver1     = list(analyte = "ergothioneine", units = "umol/L", specimen = "tissue", verified = TRUE),
    int_liver1    = list(analyte = "ergothioneine", units = "umol/L", specimen = "tissue", verified = TRUE),
    is_liver2     = list(analyte = "ergothioneine", units = "umol/L", specimen = "tissue", verified = TRUE),
    int_liver2    = list(analyte = "ergothioneine", units = "umol/L", specimen = "tissue", verified = TRUE),
    is_liver3     = list(analyte = "ergothioneine", units = "umol/L", specimen = "tissue", verified = TRUE),
    int_liver3    = list(analyte = "ergothioneine", units = "umol/L", specimen = "tissue", verified = TRUE),
    is_liver4     = list(analyte = "ergothioneine", units = "umol/L", specimen = "tissue", verified = TRUE),
    int_liver4    = list(analyte = "ergothioneine", units = "umol/L", specimen = "tissue", verified = TRUE),
    is_liver5     = list(analyte = "ergothioneine", units = "umol/L", specimen = "tissue", verified = TRUE),
    int_liver5    = list(analyte = "ergothioneine", units = "umol/L", specimen = "tissue", verified = TRUE),
    is_skin       = list(analyte = "ergothioneine", units = "umol/L", specimen = "tissue", verified = TRUE),
    int_skin      = list(analyte = "ergothioneine", units = "umol/L", specimen = "tissue", verified = TRUE),
    is_muscle     = list(analyte = "ergothioneine", units = "umol/L", specimen = "tissue", verified = TRUE),
    int_muscle    = list(analyte = "ergothioneine", units = "umol/L", specimen = "tissue", verified = TRUE),
    is_adipose    = list(analyte = "ergothioneine", units = "umol/L", specimen = "tissue", verified = TRUE),
    int_adipose   = list(analyte = "ergothioneine", units = "umol/L", specimen = "tissue", verified = TRUE),
    rbc_precursor = list(analyte = "ergothioneine", units = "umol", specimen = "blood cell", verified = TRUE),
    rbc_egt       = list(analyte = "ergothioneine", units = "umol/L", specimen = "blood cell", verified = TRUE),
    renal_duct    = list(analyte = "ergothioneine", units = "umol/L", specimen = "urine", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 128,
    n_studies = 2,
    age_range = "30-59 years (Study 1); 40-64 years (Study 2)",
    age_mean = "42.3 years (Study 1); 53.4 years (Study 2)",
    weight_mean = "63.5 kg (Study 1); 61.5 kg (Study 2)",
    sex_female_pct = 46,
    race_ethnicity = c(Asian = 100),
    disease_state = "healthy volunteers",
    dose_range = "5, 8, 10 and 20 mg/day oral EGT once daily for 8-16 weeks",
    regions = "Japan",
    notes = paste(
      "Study 1 (UMIN000041910) randomised 28 healthy Japanese employees",
      "aged 30-59 to 5 or 10 mg/day for 8 weeks (14 per arm; 8/28 female;",
      "baseline Supplementary Table S3). Study 2 (UMIN000044580)",
      "randomised 100 healthy Japanese adults aged 40-64 who reported",
      "memory or attention decline to 8 mg/day or placebo for 16 weeks",
      "(50 per arm, one placebo withdrawal; 51/100 female; baseline",
      "Table 3); only the 50 EGT recipients contribute PK. The 20 mg/day",
      "plasma and red-cell profiles used in the fit come from the earlier",
      "Katsube 2022 study and are not part of the 128 subjects counted",
      "here. n_subjects is the 28 Study 1 plus 100 Study 2 participants.",
      "All physiological constants are for a typical 70 kg adult",
      "(Supplementary Table S5 footnote a), which is above the observed",
      "mean weight in both studies."
    )
  )

  ini({
    # -----------------------------------------------------------------
    # Estimated parameters -- Supplementary Table S5, the paper's final
    # ("blush-up") fit of plasma and red-cell profiles at 5, 8, 10 and
    # 20 mg/day. Values are +/- SD; CV was 8%-20% (Results 3.2.5).
    # The earlier Table 1 fit, which excluded the 8 mg/day arm and
    # generated the Table 2 dose-selection simulation, gave R = 1.36,
    # VmaxD = 121 and VmaxP = 11.1; the vignette reproduces Table 2 by
    # overriding these three back to those values.
    # -----------------------------------------------------------------
    lrdiet <- log(1.52)
    label("Daily dietary ergothioneine ingestion rate (umol/h)")                  # Table S5: R 1.52 +/- 0.31 umol/h
    lvmax_duct <- log(105)
    label("Maximum rate of renal tubular reabsorption (umol/h)")                  # Table S5: VmaxD 105 +/- 14 umol/h
    lvmax_precursor <- log(11.0)
    label("Maximum rate of uptake into the red-cell precursor pool (umol/h)")     # Table S5: VmaxP 11.0 +/- 0.9 umol/h

    # -----------------------------------------------------------------
    # Fixed transporter constants -- Table 1.
    # -----------------------------------------------------------------
    lka <- fixed(log(6.00))
    label("Gastric emptying / absorption rate constant (1/h)")                    # Table 1: kab 6.00 /h, Fixed
    lkm <- fixed(log(21.0))
    label("Michaelis constant for OCTN1-mediated transport (umol/L)")             # Table 1: Km 21.0 uM, Fixed (Grundemann 2005)
    lvmax_liver <- fixed(log(10227))
    label("Maximum rate of hepatic OCTN1-mediated uptake (umol/h)")               # Table 1: VmaxH 10,227 umol/h, Fixed

    # -----------------------------------------------------------------
    # Fixed volumes (L) -- Table 1, all from Davies & Morris 1993 except
    # the renal duct (Sarashina 2020). Every tissue's extracellular space
    # is exactly 20% of its intracellular space in this table.
    # -----------------------------------------------------------------
    lvc <- fixed(log(3.00))
    label("Plasma volume (L)")                                                    # Table 1: Vp 3.00 L, Fixed
    lv_blood <- fixed(log(5.20))
    label("Whole-blood volume, used only to derive haematocrit (L)")              # Table 1: Blood 5.20 L, Fixed
    lv_rbc <- fixed(log(2.29))
    label("Red-blood-cell volume (L)")                                            # Table 1: VRBC 2.29 L, Fixed
    lv_liver <- fixed(log(1.69))
    label("Liver intracellular volume (L)")                                       # Table 1: Vh 1.69 L, Fixed
    lv_liver_ex <- fixed(log(0.338))
    label("Liver extracellular volume (L)")                                       # Table 1: Vhe 0.338 L, Fixed
    lv_muscle <- fixed(log(35.0))
    label("Muscle intracellular volume (L)")                                      # Table 1: Vm 35.0 L, Fixed
    lv_muscle_ex <- fixed(log(7.00))
    label("Muscle extracellular volume (L)")                                      # Table 1: Vme 7.00 L, Fixed
    lv_skin <- fixed(log(7.8))
    label("Skin intracellular volume (L)")                                        # Table 1: Vs 7.8 L, Fixed
    lv_skin_ex <- fixed(log(1.56))
    label("Skin extracellular volume (L)")                                        # Table 1: Vse 1.56 L, Fixed
    lv_adipose <- fixed(log(10.0))
    label("Adipose intracellular volume (L)")                                     # Table 1: Va 10.0 L, Fixed
    lv_adipose_ex <- fixed(log(2.00))
    label("Adipose extracellular volume (L)")                                     # Table 1: Vae 2.00 L, Fixed
    lv_duct <- fixed(log(0.127))
    label("Renal tubular duct volume (L)")                                        # Table 1: Vd 0.127 L, Fixed (Sarashina 2020)

    # -----------------------------------------------------------------
    # Fixed plasma flows (L/h) -- Table 1 (Davies & Morris 1993).
    # -----------------------------------------------------------------
    lq_liver <- fixed(log(48.7))
    label("Hepatic plasma flow (L/h)")                                            # Table 1: Qh 48.7 L/h, Fixed
    lq_muscle <- fixed(log(45.0))
    label("Muscle plasma flow (L/h)")                                             # Table 1: Qm 45.0 L/h, Fixed
    lq_skin <- fixed(log(18.0))
    label("Skin plasma flow (L/h)")                                               # Table 1: Qs 18.0 L/h, Fixed
    lq_adipose <- fixed(log(8.74))
    label("Adipose plasma flow (L/h)")                                            # Table 1: Qa 8.74 L/h, Fixed
    lq_urine <- fixed(log(0.0583))
    label("Urinary flow rate (L/h)")                                              # Table 1: Qu 0.0583 L/h, Fixed
    lgfr <- fixed(log(7.50))
    label("Glomerular filtration rate (L/h)")                                     # Table 1: GFR 7.50 L/h, Fixed

    # -----------------------------------------------------------------
    # Fixed red-cell turnover -- Table 1.
    # -----------------------------------------------------------------
    lk_rbc_diff <- fixed(log(0.0104))
    label("Red-cell differentiation rate constant (1/h)")                         # Table 1: k1 0.0104 /h, Fixed (Yasukouchi 1999)
    lk_rbc_deg <- fixed(log(0.000963))
    label("Red-cell degradation rate constant (1/h)")                             # Table 1: kb 0.000963 /h, Fixed (Obara 1962)

    # -----------------------------------------------------------------
    # Fixed tissue-to-plasma concentration ratios -- Table 1, measured in
    # mice (Araragi et al., manuscript in preparation) and restated in
    # Methods 2.1.3. These set each tissue's efflux clearance via
    # PSeff = Rt * VmaxH / Km / Kp.
    # -----------------------------------------------------------------
    lkp_liver <- fixed(log(80.1))
    label("Liver-to-plasma concentration ratio (unitless)")                       # Table 1: KpH 80.1, Fixed
    lkp_muscle <- fixed(log(5.06))
    label("Muscle-to-plasma concentration ratio (unitless)")                      # Table 1: KpM 5.06, Fixed
    lkp_skin <- fixed(log(20.8))
    label("Skin-to-plasma concentration ratio (unitless)")                        # Table 1: KpS 20.8, Fixed
    lkp_adipose <- fixed(log(17.5))
    label("Adipose-to-plasma concentration ratio (unitless)")                     # Table 1: KpA 17.5, Fixed

    # -----------------------------------------------------------------
    # Tissue-to-liver Vmax ratios (Rt in the supplement's nomenclature).
    # NOT TABULATED ANYWHERE IN THE PAPER -- derived here from the rule
    # the paper states in Methods 2.1.3: per-gram Vmax is proportional to
    # the mean OCTN1 mRNA level across three expression databases (liver
    # 0.147, muscle 0.464, adipose 0.209, skin 0.131), so the whole-organ
    # ratio is (mRNA_tissue * mass_tissue) / (mRNA_liver * mass_liver)
    # with organ mass taken from the Table 1 volumes. That choice is
    # basis-free: extracellular volume is exactly 20% of intracellular
    # volume for all four tissues, so using total or intracellular volume
    # gives identical ratios. The vignette shows the model output is
    # numerically insensitive to this derivation (a 20-fold change in
    # lr_muscle moves the Table 2 predictions by <0.01%), because uptake
    # and efflux share the same Rt factor and so the steady-state
    # tissue:plasma ratio is fixed by Kp alone.
    # -----------------------------------------------------------------
    lr_muscle <- fixed(log(65.3705))
    label("Muscle-to-liver maximum-uptake-rate ratio (unitless)")                 # derived: (0.464*42.0)/(0.147*2.028), Methods 2.1.3 + Table 1
    lr_skin <- fixed(log(4.11303))
    label("Skin-to-liver maximum-uptake-rate ratio (unitless)")                   # derived: (0.131*9.36)/(0.147*2.028), Methods 2.1.3 + Table 1
    lr_adipose <- fixed(log(8.41283))
    label("Adipose-to-liver maximum-uptake-rate ratio (unitless)")                # derived: (0.209*12.0)/(0.147*2.028), Methods 2.1.3 + Table 1

    # -----------------------------------------------------------------
    # Baseline (pre-supplementation) concentrations. These are inputs to
    # the initial-condition formulae, not estimates: the paper sets them
    # to the observed blank value of whichever arm is being simulated.
    # The defaults are the 8 mg/day row of Supplementary Table S5, which
    # is the arm the final parameter estimates above were fitted with.
    # Override for other arms -- Table 1 lists 3.42/624 (5 mg/day),
    # 3.51/562 (10 mg/day) and 2.97/594 (20 mg/day), and the Table 2
    # dose-selection simulation uses 3.16 with 594. Methods 2.1.3 calls
    # 3.16 "the average of the blank values obtained at 5, 10, and
    # 20 mg/day", but the simple mean of those three is 3.30; see the
    # vignette Errata.
    # -----------------------------------------------------------------
    cp0 <- fixed(4.07)
    label("Baseline plasma ergothioneine concentration (umol/L)")                 # Table S5, 8 mg/day: CPlasma 4.07 uM
    crbc0 <- fixed(594)
    label("Baseline red-cell ergothioneine concentration (umol/L)")               # Table S5, 8 mg/day: CRBC 594 uM

    # -----------------------------------------------------------------
    # Napp fits by nonlinear least squares on mean profiles: it estimates
    # neither between-subject variability nor a residual-error variance,
    # and none is reported anywhere in Okumura 2025 or its supplement.
    # nlmixr2 model definitions require a residual-error term, so the two
    # below are fixed placeholders for syntactic completeness only and
    # must NOT be read as estimates. Same convention as
    # Aoki_2024_bosentan_pbpk, Mi_2023_cefquinome_pbpk and
    # An_2012_mitoxantrone_human_pbpk.
    # -----------------------------------------------------------------
    propSd <- fixed(0.10)
    label("Proportional residual error placeholder, plasma (fraction)")           # not reported in Okumura 2025; placeholder only
    propSd_Crbc <- fixed(0.10)
    label("Proportional residual error placeholder, red cells (fraction)")        # not reported in Okumura 2025; placeholder only
  })

  model({
    # 1. Back-transform to the natural scale.
    rdiet <- exp(lrdiet)
    vmax_duct <- exp(lvmax_duct)
    vmax_precursor <- exp(lvmax_precursor)
    ka <- exp(lka)
    km <- exp(lkm)
    vmax_liver <- exp(lvmax_liver)
    vc <- exp(lvc)
    v_blood <- exp(lv_blood)
    v_rbc <- exp(lv_rbc)
    v_liver <- exp(lv_liver)
    v_liver_ex <- exp(lv_liver_ex)
    v_muscle <- exp(lv_muscle)
    v_muscle_ex <- exp(lv_muscle_ex)
    v_skin <- exp(lv_skin)
    v_skin_ex <- exp(lv_skin_ex)
    v_adipose <- exp(lv_adipose)
    v_adipose_ex <- exp(lv_adipose_ex)
    v_duct <- exp(lv_duct)
    q_liver <- exp(lq_liver)
    q_muscle <- exp(lq_muscle)
    q_skin <- exp(lq_skin)
    q_adipose <- exp(lq_adipose)
    q_urine <- exp(lq_urine)
    gfr <- exp(lgfr)
    k_rbc_diff <- exp(lk_rbc_diff)
    k_rbc_deg <- exp(lk_rbc_deg)
    kp_liver <- exp(lkp_liver)
    kp_muscle <- exp(lkp_muscle)
    kp_skin <- exp(lkp_skin)
    kp_adipose <- exp(lkp_adipose)
    r_muscle <- exp(lr_muscle)
    r_skin <- exp(lr_skin)
    r_adipose <- exp(lr_adipose)

    # 2. Derived terms.
    #    The liver is five equal tandem sub-compartments, so each holds a
    #    fraction fdisp = 1/5 of the hepatic volume and receives fdisp of
    #    the liver-wide uptake and efflux clearances. The Napp code writes
    #    this as an explicit `/5` on VmaxH and PSeffH in every hepatic ODE.
    fdisp <- 0.2

    #    Whole-organ maximum uptake rates, scaled off the liver by the
    #    OCTN1 expression ratios (Supporting Information S1, "Preparative
    #    calculation": VmaxS = Rs*VmaxH, VmaxM = Rm*VmaxH, VmaxA = Ra*VmaxH).
    vmax_muscle <- r_muscle * vmax_liver
    vmax_skin <- r_skin * vmax_liver
    vmax_adipose <- r_adipose * vmax_liver

    #    Intrinsic efflux clearances (L/h). Same block: PSeffT is set so
    #    that the linear-regime uptake:efflux ratio (Vmax/Km)/PSeff equals
    #    the tissue-to-plasma ratio Kp.
    pseff_liver <- vmax_liver / km / kp_liver
    pseff_muscle <- r_muscle * vmax_liver / km / kp_muscle
    pseff_skin <- r_skin * vmax_liver / km / kp_skin
    pseff_adipose <- r_adipose * vmax_liver / km / kp_adipose

    #    Haematocrit, used only to report whole-blood EGT via the paper's
    #    Equation 1, Cblood = Hct*Crbc + (1-Hct)*Cplasma.
    hct <- v_rbc / v_blood

    #    Steady-state hepatic initial concentrations (Supporting
    #    Information S1, "Estimation of initial values").
    che0 <- rdiet / q_liver + cp0
    chi0 <- vmax_liver * che0 / (km + che0) / pseff_liver

    # 3. ODE system -- transcribed from the Napp code block of Supporting
    #    Information S1. The tissue states hold CONCENTRATIONS (umol/L);
    #    `depot` and `rbc_precursor` hold AMOUNTS (umol).

    #    Gastrointestinal tract: continuous dietary intake plus the
    #    supplement dose, emptying first-order into the portal inflow.
    d/dt(depot) <- rdiet - ka * depot

    #    Plasma: return from the last hepatic sub-compartment and from the
    #    three peripheral extracellular spaces, saturable loss into the
    #    red-cell precursor pool, return from degrading red cells,
    #    glomerular filtration out, and saturable tubular reabsorption in.
    d/dt(central) <- (-(q_liver + q_muscle + q_skin + q_adipose) * central +
      q_liver * is_liver5 + q_muscle * is_muscle + q_skin * is_skin +
      q_adipose * is_adipose -
      vmax_precursor * central / (km + central) +
      k_rbc_deg * rbc_egt * v_rbc -
      gfr * central +
      vmax_duct * renal_duct / (km + renal_duct)) / vc

    #    Hepatic tandem chain. Sub-compartment 1 additionally receives the
    #    absorbed dose (complete gastrointestinal absorption is assumed);
    #    each subsequent one receives its predecessor's outflow, so
    #    is_liver5 is what returns to plasma above.
    d/dt(is_liver1) <- (ka * depot + q_liver * (central - is_liver1) -
      fdisp * vmax_liver * is_liver1 / (km + is_liver1) +
      fdisp * pseff_liver * int_liver1) / (fdisp * v_liver_ex)
    d/dt(int_liver1) <- (fdisp * vmax_liver * is_liver1 / (km + is_liver1) -
      fdisp * pseff_liver * int_liver1) / (fdisp * v_liver)

    d/dt(is_liver2) <- (q_liver * (is_liver1 - is_liver2) -
      fdisp * vmax_liver * is_liver2 / (km + is_liver2) +
      fdisp * pseff_liver * int_liver2) / (fdisp * v_liver_ex)
    d/dt(int_liver2) <- (fdisp * vmax_liver * is_liver2 / (km + is_liver2) -
      fdisp * pseff_liver * int_liver2) / (fdisp * v_liver)

    d/dt(is_liver3) <- (q_liver * (is_liver2 - is_liver3) -
      fdisp * vmax_liver * is_liver3 / (km + is_liver3) +
      fdisp * pseff_liver * int_liver3) / (fdisp * v_liver_ex)
    d/dt(int_liver3) <- (fdisp * vmax_liver * is_liver3 / (km + is_liver3) -
      fdisp * pseff_liver * int_liver3) / (fdisp * v_liver)

    d/dt(is_liver4) <- (q_liver * (is_liver3 - is_liver4) -
      fdisp * vmax_liver * is_liver4 / (km + is_liver4) +
      fdisp * pseff_liver * int_liver4) / (fdisp * v_liver_ex)
    d/dt(int_liver4) <- (fdisp * vmax_liver * is_liver4 / (km + is_liver4) -
      fdisp * pseff_liver * int_liver4) / (fdisp * v_liver)

    d/dt(is_liver5) <- (q_liver * (is_liver4 - is_liver5) -
      fdisp * vmax_liver * is_liver5 / (km + is_liver5) +
      fdisp * pseff_liver * int_liver5) / (fdisp * v_liver_ex)
    d/dt(int_liver5) <- (fdisp * vmax_liver * is_liver5 / (km + is_liver5) -
      fdisp * pseff_liver * int_liver5) / (fdisp * v_liver)

    #    Permeability-limited peripheral tissues: perfusion exchange with
    #    plasma into the extracellular space, then saturable OCTN1 uptake
    #    into and passive efflux out of the intracellular space.
    d/dt(is_muscle) <- (q_muscle * (central - is_muscle) -
      vmax_muscle * is_muscle / (km + is_muscle) +
      pseff_muscle * int_muscle) / v_muscle_ex
    d/dt(int_muscle) <- (vmax_muscle * is_muscle / (km + is_muscle) -
      pseff_muscle * int_muscle) / v_muscle

    d/dt(is_skin) <- (q_skin * (central - is_skin) -
      vmax_skin * is_skin / (km + is_skin) +
      pseff_skin * int_skin) / v_skin_ex
    d/dt(int_skin) <- (vmax_skin * is_skin / (km + is_skin) -
      pseff_skin * int_skin) / v_skin

    d/dt(is_adipose) <- (q_adipose * (central - is_adipose) -
      vmax_adipose * is_adipose / (km + is_adipose) +
      pseff_adipose * int_adipose) / v_adipose_ex
    d/dt(int_adipose) <- (vmax_adipose * is_adipose / (km + is_adipose) -
      pseff_adipose * int_adipose) / v_adipose

    #    Erythroid lineage: saturable uptake from plasma into the marrow
    #    precursor pool (an amount, umol), first-order maturation into
    #    circulating red cells, first-order red-cell degradation returning
    #    EGT to plasma. This chain is why whole-blood EGT keeps rising for
    #    weeks after plasma has levelled off.
    d/dt(rbc_precursor) <- vmax_precursor * central / (km + central) -
      k_rbc_diff * rbc_precursor
    d/dt(rbc_egt) <- (k_rbc_diff * rbc_precursor -
      k_rbc_deg * rbc_egt * v_rbc) / v_rbc

    #    Proximal tubular lumen: glomerular filtrate in, saturable OCTN1
    #    reabsorption back to plasma, urinary loss out. This is the only
    #    elimination pathway in the model.
    d/dt(renal_duct) <- (gfr * central -
      vmax_duct * renal_duct / (km + renal_duct) -
      q_urine * renal_duct) / v_duct

    # 4. Initial conditions -- Supporting Information S1, "Estimation of
    #    initial values". The gastrointestinal, hepatic, peripheral,
    #    precursor and tubular states are placed at the steady state
    #    implied by the dietary intake rate and the baseline plasma
    #    concentration; plasma and red cells start at the observed blank
    #    values of the arm being simulated.
    depot(0) <- rdiet / ka
    central(0) <- cp0
    is_liver1(0) <- che0
    is_liver2(0) <- che0
    is_liver3(0) <- che0
    is_liver4(0) <- che0
    is_liver5(0) <- che0
    int_liver1(0) <- chi0
    int_liver2(0) <- chi0
    int_liver3(0) <- chi0
    int_liver4(0) <- chi0
    int_liver5(0) <- chi0
    is_muscle(0) <- cp0
    int_muscle(0) <- vmax_muscle * cp0 / (km + cp0) / pseff_muscle
    is_skin(0) <- cp0
    int_skin(0) <- vmax_skin * cp0 / (km + cp0) / pseff_skin
    is_adipose(0) <- cp0
    int_adipose(0) <- vmax_adipose * cp0 / (km + cp0) / pseff_adipose
    rbc_precursor(0) <- vmax_precursor * cp0 / (km + cp0) / k_rbc_diff
    rbc_egt(0) <- crbc0
    renal_duct(0) <- rdiet / q_urine

    # 5. Observations. Both tissue states already hold concentrations in
    #    umol/L. Cblood is the whole-blood concentration the assay
    #    actually measures, recovered from Equation 1 of the paper; it is
    #    reported as a derived output rather than a fitted endpoint
    #    because the fit was performed on plasma and back-calculated
    #    red-cell concentrations.
    Cc <- central
    Crbc <- rbc_egt
    Cblood <- hct * rbc_egt + (1 - hct) * central

    Cc ~ prop(propSd)
    Crbc ~ prop(propSd_Crbc)
  })
}
