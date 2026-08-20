Aoki_2024_intratarget_microdosing_pbpk <- function() {
  description <- paste(
    "PBPK-PKRO (simplified whole-body, permeability-limited target tissue).",
    "Drug-agnostic small-molecule model used by Aoki, Rowland & Sugiyama",
    "(2024) to ask when Intra-Target Microdosing (ITM) -- direct delivery of",
    "a Phase-0 microdose into the target tissue -- can engage a target as",
    "well as a systemically administered therapeutic dose. Plasma exchanges",
    "with three perfusion-limited well-stirred tissues (muscle, skin,",
    "adipose), with a permeability-limited liver (hepatic extracellular",
    "space in series with the hepatocyte), and with a permeability-limited",
    "target tissue (extracellular space in series with the intracellular",
    "space). Both the liver and the target tissue take drug up by saturable",
    "active transport plus passive diffusion, return it by passive diffusion,",
    "and metabolise it by a saturable intracellular route; renal clearance",
    "removes drug from plasma. Inside the target tissue the unbound drug",
    "binds a finite receptor pool of total amount rtot, so fractional",
    "receptor occupancy (RO) -- the paper's measure of target engagement --",
    "is an explicit model output. The same structure serves both routes: an",
    "intravenous dose enters central, an ITM dose enters the target-tissue",
    "extracellular space. This is a Monte Carlo simulation study, not a fit:",
    "the paper draws 10,000 'virtual compounds' from the log-normal",
    "parameter distributions of Table 2 and scores each one for ITM success",
    "(24-h average RO >= 60%). The values below are therefore the published",
    "MEDIANS of those distributions (Supplementary Table S1, Pctl. 50), i.e.",
    "a reference typical virtual compound rather than any single compound;",
    "every parameter is fixed() because none was estimated from data. The",
    "vignette reproduces the nine fully specified example compounds of",
    "Supplementary Table S2 and their published outcomes. See the vignette",
    "Errata for three published-table discrepancies, each settled by the",
    "authors' own code."
  )
  reference <- paste(
    "Aoki Y, Rowland M, Sugiyama Y. When to consider intra-target",
    "microdosing: physiologically based pharmacokinetic modeling approach to",
    "quantitatively identify key factors for observing target engagement.",
    "Front Pharmacol. 2024;15:1366160. doi:10.3389/fphar.2024.1366160.",
    "The ODE system is transcribed line for line from the `model_text` block",
    "of the authors' own simulation script, Supplementary Data Sheet 1",
    "(DataSheet1.ZIP -> Aoki_etal_ITM_simulationCode/suppMaterial_modelCode.R);",
    "the physiological constants are Table 1; the parameter values are the",
    "Pctl. 50 column of Supplementary Table S1 (Supplementary Data Sheet 2,",
    "DataSheet2.docx); and the nine worked example compounds are",
    "Supplementary Table S2, each of which reproduces exactly a numbered row",
    "of the companion parameter file",
    "Aoki_etal_ITM_simulationCode/SuppMaterial_modelParameter.csv.",
    "The structure originates with Koyama S, Toshimoto K, Lee W, Aoki Y,",
    "Sugiyama Y. Revisiting nonlinear bosentan pharmacokinetics by",
    "physiologically based pharmacokinetic modeling. Drug Metab Dispos.",
    "2021;49(4):298-304. doi:10.1124/dmd.120.000023 (not open access; not on",
    "disk -- but it is not needed here, because Aoki 2024 publishes the",
    "complete ODE system and every parameter value in its own supplement).",
    sep = " "
  )
  vignette <- "Aoki_2024_intratarget_microdosing_pbpk"

  # The three AUC states are plain integrators that the published code
  # carries as ODE states so that the paper's endpoints can be read
  # straight off the solve: `auc_occupancy` is the integral of fractional
  # receptor occupancy, so 24-h average RO -- the ITM-success criterion --
  # is auc_occupancy(24) / 24 exactly, rather than a trapezoidal
  # approximation on whatever observation grid the caller happens to use.
  # `auc_central` and `auc_int_tumor` are the corresponding plasma and
  # target-tissue intracellular exposure integrals (the published code's
  # AUC_RDcomplex / AUC_Central / AUC_tumourIC). They are numerical
  # conveniences specific to this paper's scoring procedure, so they are
  # declared here rather than registered as canonical compartments --
  # the same treatment `occupancy` gets in the sibling model
  # Aoki_2024_bosentan_pbpk.R.
  paper_specific_compartments <- c(
    "auc_occupancy", "auc_central", "auc_int_tumor"
  )

  # Time in hours. The published code works in molar units throughout: the
  # tissue states hold concentrations in umol/L (micromolar) and the
  # receptor states hold amounts in umol, so doses must be supplied in
  # umol. The paper assumes a molar mass of 400 g/mol for the conversion
  # between its microgram dose caps and umol (Table 2, "Molar mass ...
  # Fixed ... 400"), so 100 ug = 100 / 400 = 0.25 umol.
  units <- list(time = "h", dosing = "umol", concentration = "umol/L")

  # No covariates. Every physiological volume and blood flow is a fixed
  # function of body weight (Table 1) and the paper assumes a single body
  # weight of 75 kg for all subjects, so weight enters as the fixed
  # parameter `wt` below rather than as a data column.
  covariateData <- list()

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against the ODE system in
  # suppMaterial_modelCode.R: the eight drug-distribution states hold
  # CONCENTRATIONS (umol/L) and the two receptor states hold AMOUNTS
  # (umol), which is why every dose has to be divided by a volume on the
  # way in (see the f() terms in model()).
  # The target tissue is parameterised as a solid tumour (Table 1, after
  # Rose et al. 2019), so its two spaces take the `tumor` specimen; the two
  # hepatic spaces take the generic `tissue` specimen because the specimen
  # vocabulary does not resolve extracellular fluid from hepatocyte.
  compartmentData <- list(
    central       = list(analyte = "drug", units = "umol/L", specimen = "plasma", verified = TRUE),
    is_liver      = list(analyte = "drug", units = "umol/L", specimen = "tissue", verified = TRUE),
    int_liver     = list(analyte = "drug", units = "umol/L", specimen = "tissue", verified = TRUE),
    muscle        = list(analyte = "drug", units = "umol/L", specimen = "tissue", verified = TRUE),
    skin          = list(analyte = "drug", units = "umol/L", specimen = "tissue", verified = TRUE),
    adipose       = list(analyte = "drug", units = "umol/L", specimen = "tissue", verified = TRUE),
    is_tumor      = list(analyte = "drug", units = "umol/L", specimen = "tumor", verified = TRUE),
    int_tumor     = list(analyte = "drug", units = "umol/L", specimen = "tumor", verified = TRUE),
    target        = list(analyte = "free receptor", units = "umol", specimen = "tumor", verified = TRUE),
    complex       = list(analyte = "drug-receptor complex", units = "umol", specimen = "tumor", verified = TRUE),
    auc_occupancy = list(analyte = "fractional receptor occupancy", units = "h", specimen = "not applicable", verified = TRUE),
    auc_central   = list(analyte = "drug", units = "umol*h/L", specimen = "plasma", verified = TRUE),
    auc_int_tumor = list(analyte = "drug", units = "umol*h/L", specimen = "tumor", verified = TRUE)
  )

  population <- list(
    species    = "human",
    n_subjects = paste(
      "No human subjects. 10,000 simulated 'virtual compounds' in the main",
      "analysis (Methods, 'Summarize the success rate of the ITM'), plus",
      "1,000 per scenario x 21 scenarios = 21,000 in the sensitivity",
      "analyses. Supplementary Data Sheet 1 ships the first 1,000 draws as",
      "SuppMaterial_modelParameter.csv."
    ),
    weight_range = "75 kg assumed for every subject (Methods; Table 1)",
    dose_range = paste(
      "Intravenous doses scanned over 0.01-900,000 umol on a 1-9 x decade",
      "grid to locate each compound's 'estimated therapeutic dose' (the",
      "lowest dose sustaining 24-h average RO >= 60%). The paired ITM dose",
      "is the lesser of 1/100th of that therapeutic dose or 100 ug",
      "(= 0.25 umol at the assumed 400 g/mol), delivered into the target",
      "tissue extracellular space."
    ),
    disease_state = paste(
      "Not a disease model. The target tissue is parameterised as a 5 g",
      "solid tumour, using the vascular (0.204) and interstitial (0.296)",
      "volume fractions and the perfusion rate of the permeability-limited",
      "tumour model of Rose et al. (2019) (Table 1)."
    ),
    notes = paste(
      "This is a methodological simulation study, so there is no observed",
      "dataset, no between-subject variability and no residual-error model.",
      "The 'population' is a distribution over COMPOUNDS, not over people:",
      "kinetic parameters are drawn log-normally with the interquartile",
      "ranges of Table 2, and the correlated block (fp = fb, CLr and the",
      "tissue partition scaling Kp_scaling) is drawn from the multivariate",
      "log-normal of Kato et al. (2003), while target-binding parameters",
      "(Kd, koff, receptor abundance) follow Dahl and Akerud (2013). The",
      "authors validate the resulting virtual compounds against 164",
      "marketed small molecules from Jansson et al. (2020) (Figure 3)."
    )
  )

  ini({
    # -------------------------------------------------------------------
    # Nothing in this model was estimated from data -- it is a Monte Carlo
    # simulation study -- so every parameter is fixed(). The values are the
    # medians (Pctl. 50) of the 10,000-draw parameter distribution in
    # Supplementary Table S1, which is the paper's own summary of the
    # Table 2 distributions. They therefore describe a reference typical
    # virtual compound; because each column is a marginal median, the set
    # is not one of the 10,000 draws. The vignette carries the nine fully
    # specified example compounds of Supplementary Table S2 instead, and
    # those are the quantitative validation targets.
    # -------------------------------------------------------------------

    # --- Target binding -------------------------------------------------
    # The published code takes Kd and k_on as its inputs and derives the
    # off-rate as k_off = Kd * k_on (suppMaterial_modelCode.R, first line
    # of model_text), so those are the two parameters carried here.
    lkd <- fixed(log(0.0014))
    label("Drug-receptor equilibrium dissociation constant (umol/L)")             # Table S1 Pctl.50: Kd 0.0014; Table 2 IQR [0.0003, 0.006] from Dahl and Akerud (2013)
    lkon <- fixed(log(739))
    label("Drug-receptor association rate constant (L/(umol*h))")                 # Table S1 Pctl.50: k_on 739 (see vignette Errata on the printed '1/M hr' unit)
    lrtot <- fixed(log(0.01))
    label("Total receptor abundance in the target tissue (umol)")                  # Table S1 Pctl.50: X_TotalR 0.01; Table 2 IQR [0.001, 0.1]

    # --- Hepatic disposition --------------------------------------------
    lvmax_uptake <- fixed(log(103))
    label("Maximum rate of active hepatic uptake (umol/h)")                        # Table S1 Pctl.50: Vmax_uptake 103; Table 2 IQR [10, 1000]
    lkm_uptake <- fixed(log(0.97))
    label("Michaelis constant for active hepatic uptake, unbound (umol/L)")        # Table S1 Pctl.50: Km_uptake 0.97; Table 2 IQR [1/3, 3]
    lps_dif <- fixed(log(30))
    label("Passive-diffusion clearance, liver (L/h)")                             # Table S1 Pctl.50: PSdif_inf 30; Table 2 IQR [3, 300]
    lvmax_met <- fixed(log(259))
    label("Maximum rate of hepatic metabolism (umol/h)")                           # Table S1 Pctl.50: Vmax_met 259; Table 2 IQR [25, 2500]
    lkm_met <- fixed(log(10))
    label("Michaelis constant for hepatic metabolism, unbound (umol/L)")           # Table S1 Pctl.50: Km_met 10; Table 2 IQR [10/3, 30]

    # --- Target-tissue disposition --------------------------------------
    # Table 2 sets the target-tissue transport and metabolism capacities to
    # the hepatic values scaled to the size of the target organ: the active
    # routes are "1/5 of hepatic metabolism normalized to the size of the
    # organ" and the passive route is "the same as the passive uptake of
    # the compound to liver normalized to the size of the organ".
    lvmax_uptake_tumor <- fixed(log(0.2))
    label("Maximum rate of active uptake into the target tissue (umol/h)")         # Table S1 Pctl.50: Vmax_tumorUptake 0.2; Table 2 IQR [0.0192, 1.92]
    lkm_uptake_tumor <- fixed(log(1))
    label("Michaelis constant for active target-tissue uptake, unbound (umol/L)")  # Table S1 Pctl.50: Km_tumourUptake 1; Table 2 IQR [1/3, 3]
    lps_dif_tumor <- fixed(log(0.12))
    label("Passive-diffusion clearance, target tissue (L/h)")                      # Table S1 Pctl.50: PSdiff_tumInflux 0.12 (= PSdiff_tumEflux); Table 2 IQR [0.0115, 1.15]
    lvmax_met_tumor <- fixed(log(0.19))
    label("Maximum rate of metabolism inside the target tissue (umol/h)")          # Table S1 Pctl.50: Vmax_tumMet 0.19; Table 2 IQR [0.0192, 1.92]
    lkm_met_tumor <- fixed(log(10))
    label("Michaelis constant for target-tissue metabolism, unbound (umol/L)")     # Table S1 Pctl.50: Km_tumMet 10; Table 2 IQR [10/3, 30]

    # --- Unbound fractions ----------------------------------------------
    # Table 2 states fp = fb (plasma and blood unbound fractions are taken
    # to be equal); the CSV that accompanies the code additionally carries
    # fh and ft as exactly equal to one another for every draw.
    fu_b <- fixed(0.26)
    label("Unbound fraction in plasma / blood (unitless)")                          # Table S1 Pctl.50: fb 0.26, correlated draw from Kato et al. (2003)
    fu_liver <- fixed(0.83)
    label("Unbound fraction in the hepatocyte (unitless)")                          # Table S1 Pctl.50: fh 0.83
    fu_tumor <- fixed(0.83)
    label("Unbound fraction in the target-tissue intracellular space (unitless)")   # Table S1 Pctl.50: ft 0.83

    # --- Systemic disposition -------------------------------------------
    lvc <- fixed(log(6))
    label("Central (plasma) volume of distribution (L)")                            # Table S1: V_central fixed at 6 for all 10,000 draws
    lcl_renal <- fixed(log(1.7))
    label("Renal clearance from plasma (L/h)")                                      # Table S1 Pctl.50: CLr 1.7, correlated draw from Kato et al. (2003)
    lkp_muscle <- fixed(log(0.052))
    label("Muscle-to-plasma partition coefficient (unitless)")                       # Table S1 Pctl.50: Kpm 0.052; Table 2: Kpm = 0.2 * Kp_scaling
    lkp_skin <- fixed(log(0.13))
    label("Skin-to-plasma partition coefficient (unitless)")                          # Table S1 Pctl.50: Kps 0.13; Table 2: Kps = 0.5 * Kp_scaling
    lkp_adipose <- fixed(log(0.026))
    label("Adipose-to-plasma partition coefficient (unitless)")                       # Table S1 Pctl.50: Kpa 0.026; Table 2: Kpa = 0.1 * Kp_scaling

    # --- Target-tissue physiology ---------------------------------------
    # Table 1, from the permeability-limited tumour model of Rose et al.
    # (2019), for a 5 g target tissue. These three are carried as their own
    # parameters (rather than derived from a tumour mass) because the
    # paper's sensitivity analyses vary the target blood flow and the
    # target volume independently -- see the vignette, which overrides
    # q_tumor across seven values to reproduce Figure 4A.
    lq_tumor <- fixed(log(0.2058))
    label("Blood flow rate to the target tissue (L/h)")                              # Table 1: Qt = 0.686 * 60 / 1000 * 5 = 0.2058
    lv_tumor <- fixed(log(0.0025))
    label("Extracellular (vascular + interstitial) volume of the target tissue (L)") # Table 1: Vt = 5 / 1000 * (0.204 + 0.296) = 0.0025
    lv_tumor_ic <- fixed(log(0.0025))
    label("Intracellular volume of the target tissue (L)")                            # Table 1: Vt_ic = 5 / 1000 * (1 - 0.204 - 0.296) = 0.0025

    # --- Body weight ----------------------------------------------------
    # Every other physiological volume and flow is a linear function of
    # body weight (Table 1), evaluated in model() so that the Table 1
    # formulae stay visible and a user can rescale the whole body by
    # overriding this one parameter.
    wt <- fixed(75)
    label("Body weight (kg)")                                                        # Methods: "assuming all subjects in the trials weigh 75 kg"; Table 1 "Weight"

    # -------------------------------------------------------------------
    # There is no observed concentration data anywhere in this paper and
    # hence no residual-error model. nlmixr2 model definitions require a
    # residual-error term, so propSd is a fixed placeholder for syntactic
    # completeness only and must NOT be read as an estimate. Same
    # convention as the sibling Aoki_2024_bosentan_pbpk.R.
    # -------------------------------------------------------------------
    propSd <- fixed(0.10)
    label("Proportional residual error placeholder (fraction)")                      # not reported by Aoki 2024; placeholder only
  })

  model({
    # 1. Back-transform to the natural scale.
    kd <- exp(lkd)
    kon <- exp(lkon)
    rtot <- exp(lrtot)
    vmax_uptake <- exp(lvmax_uptake)
    km_uptake <- exp(lkm_uptake)
    ps_dif <- exp(lps_dif)
    vmax_met <- exp(lvmax_met)
    km_met <- exp(lkm_met)
    vmax_uptake_tumor <- exp(lvmax_uptake_tumor)
    km_uptake_tumor <- exp(lkm_uptake_tumor)
    ps_dif_tumor <- exp(lps_dif_tumor)
    vmax_met_tumor <- exp(lvmax_met_tumor)
    km_met_tumor <- exp(lkm_met_tumor)
    vc <- exp(lvc)
    cl_renal <- exp(lcl_renal)
    kp_muscle <- exp(lkp_muscle)
    kp_skin <- exp(lkp_skin)
    kp_adipose <- exp(lkp_adipose)
    q_tumor <- exp(lq_tumor)
    v_tumor <- exp(lv_tumor)
    v_tumor_ic <- exp(lv_tumor_ic)

    # 2. Receptor off-rate, from the published identity k_off = Kd * k_on
    #    (suppMaterial_modelCode.R, `k_off=Kd*k_on`).
    koff <- kd * kon

    # 3. Weight-scaled physiology, transcribed from Table 1 exactly as the
    #    published code writes it. Flows are given per kg per minute and
    #    converted to L/h by the * 60 / 1000 factor; volumes are given in
    #    mL/kg and converted to L by the / 1000 factor. NOTE that Table 1
    #    mislabels Va as "Volume of skin"; it is the adipose volume, as
    #    the code's pairing of Va with Kpa confirms (see vignette Errata).
    q_adipose <- 3.72 * wt * 60 / 1000
    q_liver <- 20.7 * wt * 60 / 1000
    q_muscle <- 10.7 * wt * 60 / 1000
    q_skin <- 4.28 * wt * 60 / 1000
    v_adipose <- 142 * wt / 1000
    v_liver <- 17.4 * wt / 1000
    v_liver_ex <- 6.69 * wt / 1000
    v_muscle <- 429 * wt / 1000
    v_skin <- 111 * wt / 1000

    # 4. ODE system -- transcribed line for line from the `model_text`
    #    block of Supplementary Data Sheet 1
    #    (Aoki_etal_ITM_simulationCode/suppMaterial_modelCode.R).
    #
    #    The published code names five flux intermediates (HepUptake,
    #    HepEflux, HepMet, Tumour_uptake, Tumour_eflux, Tumour_met) and
    #    reuses each in two ODEs. They are written out INLINE here, twice
    #    where necessary, on purpose: in the nlmixr2 model-function form a
    #    named intermediate that references an ODE state and is then used
    #    inside d/dt() can silently evaluate to zero, which would delete a
    #    Michaelis-Menten term without any error or warning. The vignette
    #    additionally solves the authors' verbatim RxODE model_text
    #    alongside this one and asserts the two agree, which is the direct
    #    check that no term was lost in transcription.

    #    Plasma: hepatic and target-tissue venous return, renal loss, and
    #    perfusion exchange with the three well-stirred tissues. The liver
    #    and target-tissue extracellular spaces exchange with plasma
    #    directly (no partition coefficient); muscle, skin and adipose
    #    return through their tissue-to-plasma partition coefficients.
    d/dt(central) <- (q_liver * (is_liver - central) -
      cl_renal * central +
      q_muscle * (muscle / kp_muscle - central) +
      q_skin * (skin / kp_skin - central) +
      q_adipose * (adipose / kp_adipose - central) +
      q_tumor * (is_tumor - central)) / vc

    #    Hepatic extracellular (sinusoidal) space: perfusion in and out,
    #    minus saturable active uptake plus passive influx into the
    #    hepatocyte, plus passive efflux back out of it.
    d/dt(is_liver) <- (q_liver * (central - is_liver) -
      (vmax_uptake * fu_b * is_liver / (km_uptake + fu_b * is_liver) +
        ps_dif * fu_b * is_liver) +
      ps_dif * fu_liver * int_liver) / v_liver_ex

    #    Hepatocyte: uptake in, passive efflux out, saturable metabolism.
    d/dt(int_liver) <- ((vmax_uptake * fu_b * is_liver / (km_uptake + fu_b * is_liver) +
      ps_dif * fu_b * is_liver) -
      ps_dif * fu_liver * int_liver -
      vmax_met * fu_liver * int_liver / (km_met + fu_liver * int_liver)) / v_liver

    #    Perfusion-limited, well-stirred peripheral tissues.
    d/dt(muscle) <- q_muscle * (central - muscle / kp_muscle) / v_muscle
    d/dt(skin) <- q_skin * (central - skin / kp_skin) / v_skin
    d/dt(adipose) <- q_adipose * (central - adipose / kp_adipose) / v_adipose

    #    Target-tissue extracellular space: perfusion in and out, minus
    #    saturable active uptake plus passive influx into the cell, plus
    #    passive efflux back out. This is also where an ITM dose lands.
    d/dt(is_tumor) <- (q_tumor * (central - is_tumor) -
      (vmax_uptake_tumor * fu_b * is_tumor / (km_uptake_tumor + fu_b * is_tumor) +
        ps_dif_tumor * fu_b * is_tumor) +
      ps_dif_tumor * fu_tumor * int_tumor) / v_tumor

    #    Target-tissue intracellular space: uptake in, passive efflux out,
    #    saturable metabolism, and the receptor binding / unbinding
    #    exchange. NOTE that the published code applies the PLASMA unbound
    #    fraction fu_b -- not the intracellular fu_tumor -- to the
    #    intracellular concentration in the binding term; that reading is
    #    transcribed verbatim and flagged in the vignette Errata.
    d/dt(int_tumor) <- ((vmax_uptake_tumor * fu_b * is_tumor / (km_uptake_tumor + fu_b * is_tumor) +
      ps_dif_tumor * fu_b * is_tumor) -
      ps_dif_tumor * fu_tumor * int_tumor -
      vmax_met_tumor * fu_tumor * int_tumor / (km_met_tumor + fu_tumor * int_tumor) +
      koff * complex -
      kon * fu_b * int_tumor * target) / v_tumor_ic

    #    Receptor binding layer, in AMOUNTS (umol): free receptor and the
    #    drug-receptor complex.
    d/dt(target) <- koff * complex - kon * fu_b * int_tumor * target
    d/dt(complex) <- kon * fu_b * int_tumor * target - koff * complex

    #    Exposure integrators. auc_occupancy integrates the fractional
    #    receptor occupancy, so the paper's ITM-success criterion is
    #    auc_occupancy / 24 >= 0.60 at t = 24 h.
    d/dt(auc_occupancy) <- complex / rtot
    d/dt(auc_central) <- central
    d/dt(auc_int_tumor) <- int_tumor

    # 5. Initial condition: the receptor pool starts fully unoccupied.
    #    suppMaterial_modelCode.R achieves this by dosing X_TotalR into
    #    the X_FreeR compartment at t = 0; an explicit initial condition is
    #    the equivalent idiomatic rxode2 form and spares the user's event
    #    table a receptor-priming dose row.
    target(0) <- rtot

    # 6. Dose scaling. The published code doses amount / volume into the
    #    two dosing compartments because both hold CONCENTRATIONS
    #    (`dose = iv_dose / V_central, dosing.to = 1` for the intravenous
    #    route and `dose = itm_dose / Vt, dosing.to = 7` for the ITM
    #    route). Expressing those divisions as bioavailability terms lets
    #    the user dose either compartment with an ordinary amount in umol.
    f(central) <- 1 / vc
    f(is_tumor) <- 1 / v_tumor

    # 7. Outputs. `central` already is the plasma concentration in umol/L,
    #    and fractional receptor occupancy is the bound fraction of the
    #    receptor pool -- the paper's measure of target engagement.
    occupancy <- complex / rtot
    Cc <- central
    Cc ~ prop(propSd)
  })
}
