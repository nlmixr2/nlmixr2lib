Aoki_2024_bosentan_pbpk <- function() {
  description <- paste(
    "PBPK-TMDD (semi-mechanistic, dispersion liver). Bosentan disposition",
    "after single intravenous doses of 10-750 mg in healthy adults,",
    "re-estimated by Aoki & Sugiyama (2024) with the Cluster Gauss-Newton",
    "method (CGNM) on the model structure of Koyama et al. (2021). The",
    "liver is resolved as a five-compartment tandem dispersion model: five",
    "hepatic extracellular (sinusoidal) sub-compartments in series with the",
    "hepatic blood flow, each exchanging with its own hepatocyte",
    "sub-compartment. Sinusoid-to-hepatocyte transport is the sum of",
    "saturable OATP-mediated uptake (Michaelis-Menten in unbound sinusoidal",
    "concentration) and passive diffusion; efflux back to the sinusoid is",
    "passive diffusion scaled by the influx/efflux ratio gamma_dif. Drug is",
    "eliminated by hepatic metabolism from the hepatocytes and by renal",
    "clearance from plasma. Muscle, skin and adipose are perfusion-limited",
    "well-stirred tissues. Superimposed on the PBPK backbone is a",
    "target-mediated drug disposition (TMDD) layer in plasma: unbound drug",
    "binds a finite endothelin-receptor pool of total amount rtot with",
    "dissociation constant kd and off-rate koff, which is the source of the",
    "low-dose nonlinearity and which permits receptor occupancy to be",
    "simulated. The eight parameters below without fixed() are the CGNM",
    "estimates of Aoki 2024 Table S2 (all eight identifiable when the 10 mg",
    "arm is included); every other value is a fixed physiological or",
    "compound constant taken from the published model code. CGNM is a",
    "fixed-effects nonlinear-least-squares method, so this model has no",
    "between-subject variability and the paper reports no residual-error",
    "model -- the propSd term is a placeholder. See the vignette Errata."
  )
  reference <- paste(
    "Aoki Y, Sugiyama Y. Cluster Gauss-Newton method for a quick",
    "approximation of profile likelihood: With application to",
    "physiologically-based pharmacokinetic models.",
    "CPT Pharmacometrics Syst Pharmacol. 2024;13(1):54-67.",
    "doi:10.1002/psp4.13055.",
    "Model structure, fixed physiological constants and the observed",
    "concentration data are transcribed from Supporting Information file",
    "PSP4-13-54-s004.r (the authors' RxODE / CGNM script for Example 3);",
    "the eight estimated parameters are the 'proposed method' column of",
    "Supporting Information Table S2 (file PSP4-13-54-s003.docx).",
    "The model structure originates with Koyama S, Toshimoto K, Lee W,",
    "Aoki Y, Sugiyama Y. Revisiting nonlinear bosentan pharmacokinetics by",
    "physiologically based pharmacokinetic modeling: target binding, albeit",
    "not a major contributor to nonlinearity, can offer prediction of",
    "target occupancy. Drug Metab Dispos. 2021;49(4):298-304.",
    "doi:10.1124/dmd.120.000023 (not open access; not on disk -- every",
    "value here is traced to the Aoki 2024 sources listed above).",
    sep = " "
  )
  vignette <- "Aoki_2024_bosentan_pbpk"

  # The liver is a five-compartment tandem dispersion model, so the
  # hepatic states carry a 1..5 index that no canonical compartment name
  # covers. The stems follow the canonical PBPK sub-compartment
  # vocabulary of inst/references/compartment-names.md
  # (`pbpkSubCompartmentRegex`): `is_liver` is the hepatic extracellular
  # (sinusoidal) space and `int_liver` is the intracellular
  # (hepatocyte) space; only the 1..5 dispersion index is
  # paper-specific. `occupancy` is the fractional endothelin-receptor
  # occupancy, carried as an explicit state exactly as the published
  # code carries it (`R_Occupancy` in PSP4-13-54-s004.r); it is
  # mathematically redundant with complex / rtot and the vignette
  # asserts that identity as a numerical check on the integration.
  paper_specific_compartments <- c(
    "is_liver1", "is_liver2", "is_liver3", "is_liver4", "is_liver5",
    "int_liver1", "int_liver2", "int_liver3", "int_liver4", "int_liver5",
    "occupancy"
  )

  # Time in hours. The published code works in molar units throughout:
  # the tissue states hold concentrations in umol/L (micromolar) and the
  # receptor states hold amounts in umol, so doses must be supplied in
  # umol (bosentan MW 551.61 g/mol, so mg * 1000 / 551.61 = umol).
  units <- list(time = "hour", dosing = "umol", concentration = "umol/L")

  # No covariates. Every physiological volume, flow and partition
  # coefficient in PSP4-13-54-s004.r is a fixed constant for a single
  # typical adult; the published code carries no body-weight or
  # demographic scaling of any kind.
  covariateData <- list()

  population <- list(
    species    = "human",
    dose_range = "10, 50, 250, 500 and 750 mg single intravenous doses (5 arms)",
    notes      = paste(
      "Aoki 2024 identifies the fitting data only as 'plasma",
      "concentration data from 10, 50, 250, 500, and 750 mg arms'",
      "(Results, Example 3). The 13 observations per arm are hard-coded",
      "in Supporting Information file PSP4-13-54-s004.r and are",
      "reproduced in the vignette. Subject counts, ages, weights and sex",
      "distribution are not reported by Aoki 2024; they belong to the",
      "underlying clinical study reported through Koyama 2021, which is",
      "not open access and is not on disk. The intravenous route is not",
      "stated in words by Aoki 2024 but is unambiguous from the code and",
      "the data: the dose enters the central compartment instantaneously",
      "at t = 0.08 h with no absorption compartment, and for the 750 mg",
      "arm dose / MW / vc = 750000 / 551.61 / 8.95 = 151.9 umol/L against",
      "an observed 151.1 umol/L at 0.167 h."
    )
  )

  ini({
    # -----------------------------------------------------------------
    # Estimated parameters -- Aoki 2024 Supporting Information Table S2,
    # 'proposed method' column (CGNM fit including the 10 mg arm, under
    # which all eight parameters are declared identifiable). Bracketed
    # values in that table are profile-likelihood interquartile ranges,
    # not standard errors, and are not encoded here.
    # -----------------------------------------------------------------
    lcl_met <- log(9.23)
    label("Hepatic intrinsic metabolic clearance, unbound (L/h)")                        # Table S2: CL_met 9.23 [7.7, 14]
    lkd <- log(0.0108)
    label("Endothelin-receptor dissociation equilibrium constant (umol/L)")              # Table S2: Kd 0.0108 [0.0026, 0.03]
    lkm <- log(0.665)
    label("Michaelis constant for hepatic OATP uptake, unbound (umol/L)")                # Table S2: Km_uptake 0.665 [0.58, 0.94]
    lps_dif <- log(3.5)
    label("Passive-diffusion influx clearance, sinusoid to hepatocyte (L/h)")            # Table S2: PSdif_inf 3.5 [1.5, 6.8]
    lvc <- log(8.95)
    label("Central (plasma) volume of distribution (L)")                                 # Table S2: Vcentral 8.95 [8.3, 10]
    lvmax <- log(962)
    label("Maximum rate of hepatic OATP-mediated uptake (umol/h)")                       # Table S2: Vmax_uptake 962 [640, 1500]
    lrtot <- log(20.1)
    label("Total endothelin-receptor amount (umol)")                                     # Table S2: X_TotalR 20.1 [9.6, 47]
    lkoff <- log(0.404)
    label("Endothelin-receptor dissociation rate constant (1/h)")                        # Table S2: k_off 0.404 [0.23, 0.84]

    # -----------------------------------------------------------------
    # Fixed physiological and compound constants -- the `parameter_text`
    # block at the head of Aoki 2024 Supporting Information file
    # PSP4-13-54-s004.r. These are held constant by the published code
    # (they are not among the eight elements of the CGNM parameter
    # vector), so they are wrapped in fixed().
    # -----------------------------------------------------------------
    lcl_renal <- fixed(log(0.144))
    label("Renal clearance from plasma (L/h)")                                           # s004.r parameter block: CLr=0.144
    lq_liver <- fixed(log(96.9))
    label("Hepatic (total) plasma flow (L/h)")                                           # s004.r parameter block: Qh=96.9
    lq_muscle <- fixed(log(50.1))
    label("Muscle plasma flow (L/h)")                                                    # s004.r parameter block: Qm=50.1
    lq_skin <- fixed(log(20.1))
    label("Skin plasma flow (L/h)")                                                      # s004.r parameter block: Qs=20.1
    lq_adipose <- fixed(log(17.4))
    label("Adipose plasma flow (L/h)")                                                   # s004.r parameter block: Qa=17.4
    lv_liver <- fixed(log(1.36))
    label("Total hepatocyte (intracellular) volume (L)")                                 # s004.r parameter block: Vh=1.36
    lv_liver_ex <- fixed(log(0.521))
    label("Total hepatic extracellular (sinusoidal) volume (L)")                         # s004.r parameter block: Vhe=0.521
    lv_muscle <- fixed(log(33.4))
    label("Muscle volume (L)")                                                           # s004.r parameter block: Vm=33.4
    lv_skin <- fixed(log(8.69))
    label("Skin volume (L)")                                                             # s004.r parameter block: Vs=8.69
    lv_adipose <- fixed(log(11.1))
    label("Adipose volume (L)")                                                          # s004.r parameter block: Va=11.1
    lk_muscle <- fixed(log(0.119))
    label("Muscle-to-plasma partition coefficient (unitless)")                           # s004.r parameter block: Kpm=0.119
    lk_skin <- fixed(log(0.483))
    label("Skin-to-plasma partition coefficient (unitless)")                             # s004.r parameter block: Kps=0.483
    lk_adipose <- fixed(log(0.121))
    label("Adipose-to-plasma partition coefficient (unitless)")                          # s004.r parameter block: Kpa=0.121
    fu_b <- fixed(0.033)
    label("Unbound fraction in plasma (unitless)")                                       # s004.r parameter block: fb=0.033
    fu_liver <- fixed(0.0696)
    label("Unbound fraction in hepatocyte (unitless)")                                   # s004.r parameter block: fh=0.0696
    gamma_dif <- fixed(0.243)
    label("Passive-diffusion influx/efflux ratio, sinusoid vs hepatocyte (unitless)")    # s004.r ODEs: the literal 0.243 in `PSdif_inf/0.243`

    # -----------------------------------------------------------------
    # CGNM is a fixed-effects nonlinear-least-squares method: it fits
    # the mean profile of each dose arm and estimates neither
    # between-subject variability nor a residual-error variance. The
    # objective in PSP4-13-54-s004.r is the sum of squared residuals of
    # log10(concentration) (`errormodel = function(y) log10(y)`), whose
    # nlmixr2 analogue is a proportional error, but no magnitude is
    # reported anywhere in the paper or its supplement. nlmixr2 model
    # definitions require a residual-error term, so propSd below is a
    # fixed placeholder for syntactic completeness only and must NOT be
    # read as an estimate. Same convention as Mi_2023_cefquinome_pbpk,
    # Kang_2023_artesunate_hamster_pbpk and An_2012_mitoxantrone_*_pbpk.
    # -----------------------------------------------------------------
    propSd <- fixed(0.10)
    label("Proportional residual error placeholder (fraction)")                          # not reported in Aoki 2024; placeholder only
  })

  model({
    # 1. Back-transform to the natural scale.
    cl_met <- exp(lcl_met)
    kd <- exp(lkd)
    km <- exp(lkm)
    ps_dif <- exp(lps_dif)
    vc <- exp(lvc)
    vmax <- exp(lvmax)
    rtot <- exp(lrtot)
    koff <- exp(lkoff)

    cl_renal <- exp(lcl_renal)
    q_liver <- exp(lq_liver)
    q_muscle <- exp(lq_muscle)
    q_skin <- exp(lq_skin)
    q_adipose <- exp(lq_adipose)
    v_liver <- exp(lv_liver)
    v_liver_ex <- exp(lv_liver_ex)
    v_muscle <- exp(lv_muscle)
    v_skin <- exp(lv_skin)
    v_adipose <- exp(lv_adipose)
    k_muscle <- exp(lk_muscle)
    k_skin <- exp(lk_skin)
    k_adipose <- exp(lk_adipose)

    # 2. Derived terms.
    #    The liver is split into five equal tandem sub-compartments, so
    #    each holds a fraction fdisp = 1/5 of the hepatic volume and
    #    receives fdisp of every liver-wide transport and elimination
    #    clearance. The published code writes this fraction as the
    #    literal 0.2 throughout the hepatic ODEs.
    fdisp <- 0.2
    #    Second-order association rate, from the TMDD identity
    #    kon = koff / kd (s004.r writes `k_off/Kd` inline in every
    #    binding term).
    kon <- koff / kd
    #    Bimolecular binding and unbinding fluxes (umol/h). Only unbound
    #    plasma drug binds the receptor.
    bind <- kon * fu_b * central * target
    unbind <- koff * complex
    #    Saturable OATP uptake plus passive influx, per hepatic
    #    sub-compartment (umol/h, before the fdisp weighting).
    upt1 <- vmax * fu_b * is_liver1 / (km + fu_b * is_liver1) + ps_dif * fu_b * is_liver1
    upt2 <- vmax * fu_b * is_liver2 / (km + fu_b * is_liver2) + ps_dif * fu_b * is_liver2
    upt3 <- vmax * fu_b * is_liver3 / (km + fu_b * is_liver3) + ps_dif * fu_b * is_liver3
    upt4 <- vmax * fu_b * is_liver4 / (km + fu_b * is_liver4) + ps_dif * fu_b * is_liver4
    upt5 <- vmax * fu_b * is_liver5 / (km + fu_b * is_liver5) + ps_dif * fu_b * is_liver5
    #    Passive efflux back to the sinusoid (umol/h, before fdisp).
    eff1 <- ps_dif / gamma_dif * fu_liver * int_liver1
    eff2 <- ps_dif / gamma_dif * fu_liver * int_liver2
    eff3 <- ps_dif / gamma_dif * fu_liver * int_liver3
    eff4 <- ps_dif / gamma_dif * fu_liver * int_liver4
    eff5 <- ps_dif / gamma_dif * fu_liver * int_liver5

    # 3. ODE system -- transcribed line for line from the `model_text`
    #    block of Aoki 2024 Supporting Information file
    #    PSP4-13-54-s004.r. Note that the tissue states hold
    #    CONCENTRATIONS (umol/L) while target / complex hold AMOUNTS
    #    (umol), exactly as the published code does.

    #    Plasma: hepatic return, hepatic first-pass out, renal loss,
    #    perfusion exchange with the three peripheral tissues, and the
    #    TMDD binding terms.
    d/dt(central) <- (q_liver * is_liver5 - q_liver * central -
      cl_renal * central -
      q_muscle * (central - muscle / k_muscle) -
      q_skin * (central - skin / k_skin) -
      q_adipose * (central - adipose / k_adipose) +
      unbind - bind) / vc

    #    Hepatic extracellular (sinusoidal) tandem chain. Sub-compartment
    #    1 receives arterial/portal inflow from plasma; each subsequent
    #    one receives the outflow of its predecessor, so is_liver5 is
    #    what returns to plasma above.
    d/dt(is_liver1) <- (q_liver * central - q_liver * is_liver1 -
      fdisp * upt1 + fdisp * eff1) / (fdisp * v_liver_ex)
    d/dt(is_liver2) <- (q_liver * is_liver1 - q_liver * is_liver2 -
      fdisp * upt2 + fdisp * eff2) / (fdisp * v_liver_ex)
    d/dt(is_liver3) <- (q_liver * is_liver2 - q_liver * is_liver3 -
      fdisp * upt3 + fdisp * eff3) / (fdisp * v_liver_ex)
    d/dt(is_liver4) <- (q_liver * is_liver3 - q_liver * is_liver4 -
      fdisp * upt4 + fdisp * eff4) / (fdisp * v_liver_ex)
    d/dt(is_liver5) <- (q_liver * is_liver4 - q_liver * is_liver5 -
      fdisp * upt5 + fdisp * eff5) / (fdisp * v_liver_ex)

    #    Hepatocyte (intracellular) chain: uptake in, passive efflux out,
    #    unbound intrinsic metabolic clearance out.
    d/dt(int_liver1) <- (fdisp * upt1 - fdisp * eff1 -
      fdisp * cl_met * fu_liver * int_liver1) / (fdisp * v_liver)
    d/dt(int_liver2) <- (fdisp * upt2 - fdisp * eff2 -
      fdisp * cl_met * fu_liver * int_liver2) / (fdisp * v_liver)
    d/dt(int_liver3) <- (fdisp * upt3 - fdisp * eff3 -
      fdisp * cl_met * fu_liver * int_liver3) / (fdisp * v_liver)
    d/dt(int_liver4) <- (fdisp * upt4 - fdisp * eff4 -
      fdisp * cl_met * fu_liver * int_liver4) / (fdisp * v_liver)
    d/dt(int_liver5) <- (fdisp * upt5 - fdisp * eff5 -
      fdisp * cl_met * fu_liver * int_liver5) / (fdisp * v_liver)

    #    Perfusion-limited, well-stirred peripheral tissues.
    d/dt(muscle) <- q_muscle * (central - muscle / k_muscle) / v_muscle
    d/dt(skin) <- q_skin * (central - skin / k_skin) / v_skin
    d/dt(adipose) <- q_adipose * (central - adipose / k_adipose) / v_adipose

    #    TMDD layer in plasma (amounts, umol).
    d/dt(target) <- unbind - bind
    d/dt(complex) <- bind - unbind
    d/dt(occupancy) <- (bind - unbind) / rtot

    # 4. Initial condition: the receptor pool starts fully unoccupied.
    #    PSP4-13-54-s004.r achieves this with a dose of X_TotalR into the
    #    X_FreeR compartment at t = 0; an explicit initial condition is
    #    the equivalent idiomatic rxode2 form and does not require the
    #    user's event table to carry a receptor-priming dose row.
    target(0) <- rtot

    # 5. Dose scaling. The published code doses the central compartment
    #    with `dose_amount * 1000 / MW_bos / Vcentral`, i.e. it converts
    #    the mg dose to umol and then pre-divides by the central volume
    #    because `central` is a CONCENTRATION state. Expressing that
    #    division as a bioavailability term lets the user dose `central`
    #    with an ordinary amount in umol.
    f(central) <- 1 / vc

    # 6. Observation. `central` already is the plasma concentration in
    #    umol/L.
    Cc <- central
    Cc ~ prop(propSd)
  })
}
