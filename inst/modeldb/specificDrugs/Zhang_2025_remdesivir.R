Zhang_2025_remdesivir <- function() {
  description <- paste(
    "Six-state mechanism-based cascade model for intravenous remdesivir",
    "and its two circulating metabolites GS-704277 (the paper's 'IM') and",
    "GS-441524 (the paper's 'NUC') in healthy adults and adults with renal",
    "impairment (Zhang 2025). Each of the three analytes occupies a central",
    "and a peripheral compartment; metabolism runs sequentially down the",
    "cascade in BOTH compartments, which is the feature the authors added to",
    "reproduce the second, slower decay phase of remdesivir and GS-704277",
    "that a plain two-compartment model could not capture. Remdesivir is",
    "eliminated from its central compartment only, and the peripheral",
    "GS-441524 pool drains to the intracellular active triphosphate",
    "GS-443902, which is not itself a state because no triphosphate",
    "concentration data were available. The system is written entirely in",
    "CONCENTRATION space: every state is a concentration in ng/mL and every",
    "parameter is a first-order rate constant in 1/h, so no volume of",
    "distribution appears anywhere in the model and a dose must be supplied",
    "as the initial central remdesivir concentration rather than as a mass.",
    "Renal impairment enters as a single binary covariate pooling every",
    "impairment stratum. Fitted by SAEM in Monolix to DIGITISED MEAN",
    "concentration profiles, not to individual patient data.")
  reference <- "Zhang S, Jeong S, Jiang B, Ho H. Pharmacokinetic simulations for remdesivir and its metabolites in healthy subjects and patients with renal impairment. Front Pharmacol. 2025;16:1488961. doi:10.3389/fphar.2025.1488961"
  vignette <- "Zhang_2025_remdesivir"

  # `dosing` is deliberately ng/mL, not mg. Zhang 2025 Equations 1-6 are
  # written on concentrations and Table 1 reports every parameter in 1/h;
  # no volume of distribution is reported anywhere in the paper, so the
  # published model cannot convert an administered mass into a
  # concentration. A dose therefore enters this model as the initial
  # central remdesivir CONCENTRATION. See the `compartmentData` note.
  units <- list(time = "h", dosing = "ng/mL", concentration = "ng/mL")

  # Every state holds a CONCENTRATION (ng/mL), not an amount. This is a
  # faithful transcription of Zhang 2025 Equations 1-6, whose dependent
  # variables the paper defines as "the concentrations of RDV and its
  # metabolites either in central (C_C,i) or in peripheral (C_P,i)
  # compartments". Because the two compartments of each analyte exchange
  # through the SAME rate constant in both directions (Q_i multiplies
  # (C_P - C_C) in one equation and (C_C - C_P) in the other), the
  # formulation implicitly assumes equal central and peripheral volumes;
  # the paper's Q_i and CL_C,RDV are therefore micro-constants (Q/V and
  # CL/V) rather than clearances, despite their names in Table 1.
  # The peripheral compartments are recorded as `blood cell` rather than the
  # generic `tissue` because Zhang 2025 identifies them explicitly with
  # peripheral blood mononuclear cells: the Introduction states that "RDV
  # enters peripheral blood mononuclear cells (PBMCs) and undergoes hydrolysis
  # by esterase to form a transient intermediate metabolite (GS-704277)", and
  # the Discussion interprets the peripheral metabolic terms as "the conversion
  # of IM to NUC in peripheral blood mononuclear cells". Only the central
  # (plasma) concentrations were fitted to data; the peripheral states are
  # unobserved.
  compartmentData <- list(
    central                = list(analyte = "remdesivir (GS-5734)", units = "ng/mL", specimen = "plasma", verified = TRUE),
    peripheral1            = list(analyte = "remdesivir (GS-5734)", units = "ng/mL", specimen = "blood cell", verified = TRUE),
    central_gs704277       = list(analyte = "GS-704277", units = "ng/mL", specimen = "plasma", verified = TRUE),
    peripheral1_gs704277   = list(analyte = "GS-704277", units = "ng/mL", specimen = "blood cell", verified = TRUE),
    central_gs441524       = list(analyte = "GS-441524", units = "ng/mL", specimen = "plasma", verified = TRUE),
    peripheral1_gs441524   = list(analyte = "GS-441524", units = "ng/mL", specimen = "blood cell", verified = TRUE)
  )

  covariateData <- list(
    RENALIMP = list(
      description        = "Renal impairment of any degree.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = matched healthy control",
      notes              = "Zhang 2025 Methods 2.3: 'Renal impairment was included as a categorical covariate in the model, with 0 representing controls and 1 indicating renal impairment. To maintain simplicity, varying degrees of renal impairment were not further distinguished as separate covariates.' The underlying Zhang 2020 phase I trial classified impairment by eGFR as mild (60-89 mL/min/1.73 m^2), moderate (30-59), severe (15-29) and kidney failure (< 15), and all four strata are pooled into RENALIMP = 1 here. The covariate values encoded in this model are the healthy-control and SEVERE-renal-impairment columns of Zhang 2025 Table 1, so RENALIMP = 1 reproduces the severe stratum specifically; the paper reports no separate parameter column for the mild or moderate strata. Zhang 2025 Results 3.1 names KC,IM, KP,IM and KP,NUC as the parameters the covariate was assigned to, but Table 1 additionally differs between the two columns in QIM, QNUC, KC,NUC and KP,NTP; Table 1 is transcribed here in full because its caption states that the simulation results were derived from those values. See the vignette Errata.",
      source_name        = "Renal Impairment"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 43L,
    n_studies      = 1L,
    age_range      = "Not reported in Zhang 2025.",
    disease_state  = "Adults with varying degrees of renal impairment and matched healthy controls, enrolled in the Gilead phase I, open-label, parallel-group renal-impairment study of Zhang et al. (2020): mild (n = 12), moderate (n = 11) or severe (n = 10) renal impairment, and kidney failure (n = 6 on dialysis, n = 4 without dialysis). Healthy matched controls served as the reference group; Zhang 2025 does not report the control group size, so n_subjects counts only the 43 renally impaired participants.",
    renal_function = "Classified by eGFR: mild 60-89 mL/min/1.73 m^2, moderate 30-59, severe 15-29, kidney failure < 15. Remdesivir is contraindicated below eGFR 30 mL/min/1.73 m^2 in labelling, so the severe and kidney-failure strata are off-label exposures studied specifically to characterise them.",
    dose_range     = "Single intravenous doses, assigned by impairment stratum (Zhang 2025 Results 3.1): 100 mg for mild or moderate impairment, 40 mg for severe impairment or predialysis kidney failure, and 20 mg for postdialysis or non-dialysis kidney failure.",
    notes          = "IMPORTANT PROVENANCE LIMITATION. Zhang 2025 did not have access to individual patient data: Methods 2.1 states that the profiles 'were digitized using the open-source software Engauge Digitizer (Version 12.1)' and that 'Due to the absence of individual patient-level data, only the mean drug concentration values from the profile curves were extracted.' The mixed-effects model was therefore fitted by SAEM (Monolix, Lixoft) to DIGITISED MEAN concentration-time curves, one per analyte per renal stratum, rather than to individual observations. Any between-subject variance the fit produced describes scatter between mean curves, not between patients, which is why the parameter estimates in Table 1 are reported without random-effect magnitudes and why no IIV is carried in this model file. Data source: Zhang S, Humeniuk R, Ling J, et al., the Gilead phase I renal-impairment study cited by Zhang 2025 as 'Zhang et al. (2020)'."
  )

  ini({
    # ------------------------------------------------------------------
    # STRUCTURAL PARAMETERS -- healthy-control column of Zhang 2025
    # Table 1 ("Parameters for healthy control (Zhang et al., 2020)").
    #
    # Every parameter in Table 1 carries units of /h. They are first-order
    # rate constants acting on concentrations, so they are named here for
    # the micro-constants they are (`lk12`, `lkel`, `lkmet_*`) rather than
    # for the clearance-flavoured symbols Table 1 uses (Q_i, CL_C,RDV).
    # Table 1 reports point estimates only -- no standard errors, no
    # relative standard errors, no confidence intervals -- so none of the
    # uncertainty normally recorded alongside an estimate is available.
    # ------------------------------------------------------------------

    lk12 <- log(0.19)
    label("Remdesivir central-peripheral exchange rate constant (1/h)")                      # Zhang 2025 Table 1 row Q_RDV, healthy-control column: 0.19 /h

    lkel <- log(2.99)
    label("Remdesivir elimination rate constant from the central compartment (1/h)")         # Zhang 2025 Table 1 row CL_C,RDV, healthy-control column: 2.99 /h. Table 1 calls this "Clearance of RDV" but its units are /h, so in this concentration-space formulation it is an elimination rate constant.

    lk12_gs704277 <- log(12.53)
    label("GS-704277 central-peripheral exchange rate constant (1/h)")                       # Zhang 2025 Table 1 row Q_IM, healthy-control column: 12.53 /h

    lk12_gs441524 <- log(0.039)
    label("GS-441524 central-peripheral exchange rate constant (1/h)")                       # Zhang 2025 Table 1 row Q_NUC, healthy-control column: 0.039 /h

    # Metabolic conversion rate constants. Each is estimated separately in
    # the central and in the peripheral compartment, which is the paper's
    # central structural claim: Methods 2.2 states that "Additional
    # metabolism was proposed to occur from the peripheral remdesivir
    # compartment to its respective peripheral compartment", and the
    # Discussion explains this was done to reproduce "the two-phase slower
    # decay observed after the initial rapid drop in concentrations for
    # both RDV and IM, as illustrated in Figure 6". The `_central` /
    # `_peripheral1` suffixes are stratum suffixes in the sense of
    # parameter-names.md: the same quantity estimated once per compartment,
    # so both carry a suffix and neither keeps the bare canonical name.

    lkmet_gs704277_central <- log(0.22)
    label("Remdesivir-to-GS-704277 conversion rate constant, central compartment (1/h)")     # Zhang 2025 Table 1 row K_C,IM, healthy-control column: 0.22 /h

    lkmet_gs704277_peripheral1 <- log(0.31)
    label("Remdesivir-to-GS-704277 conversion rate constant, peripheral compartment (1/h)")  # Zhang 2025 Table 1 row K_P,IM, healthy-control column: 0.31 /h

    lkmet_gs441524_central <- log(0.38)
    label("GS-704277-to-GS-441524 conversion rate constant, central compartment (1/h)")      # Zhang 2025 Table 1 row K_C,NUC, healthy-control column: 0.38 /h

    lkmet_gs441524_peripheral1 <- log(2.44)
    label("GS-704277-to-GS-441524 conversion rate constant, peripheral compartment (1/h)")   # Zhang 2025 Table 1 row K_P,NUC, healthy-control column: 2.44 /h

    # GS-441524 to the active triphosphate GS-443902, in the peripheral
    # compartment only. The estimate is enormous (3.28e10 /h) and is
    # transcribed exactly as printed. It is NOT a typo for 3.28e-10: at
    # 3.28e-10 the peripheral GS-441524 pool is never drained, the central
    # GS-441524 profile rises monotonically for a week and never peaks, and
    # Zhang 2025 Figure 3c -- which shows a clear GS-441524 peak followed
    # by a decline -- cannot be reproduced. At the printed value the
    # peripheral GS-441524 pool is driven to essentially zero
    # instantaneously, so the term acts as an irreversible sink and the
    # published figures are reproduced. The resulting system is formally
    # stiff, but LSODA handles it without difficulty: peak heights and
    # times agree to every printed digit between rxode2's default
    # tolerances and atol = 1e-12 / rtol = 1e-10, so no special solver
    # settings are required (verified in the vignette).
    lkmet_gs443902_peripheral1 <- log(3.28e10)
    label("GS-441524-to-GS-443902 conversion rate constant, peripheral compartment (1/h)")   # Zhang 2025 Table 1 row K_P,NTP, healthy-control column: 3.28E10 /h

    # ------------------------------------------------------------------
    # RENAL-IMPAIRMENT COVARIATE EFFECTS
    #
    # Zhang 2025 does not print covariate coefficients. It prints two
    # complete parameter columns in Table 1 -- one for healthy controls and
    # one for severe renal impairment -- and states in the Table 1 caption
    # that "The simulation results derived from these values are presented
    # in the Results section". Each coefficient below is therefore written
    # as the log of the ratio of the two printed column entries, which
    # keeps BOTH source values visible at the point of use and makes the
    # model reproduce Table 1 exactly at RENALIMP = 0 and RENALIMP = 1.
    #
    # Q_RDV (0.19 -> 0.19) and CL_C,RDV (2.99 -> 2.99) are identical in the
    # two columns, so they carry no covariate term at all.
    #
    # Zhang 2025 Results 3.1 and the Discussion name only KC,IM, KP,IM and
    # KP,NUC as covariate-carrying. Table 1 nonetheless differs in QIM,
    # QNUC, KC,NUC and KP,NTP as well. Table 1 is transcribed in full here
    # because it is the numerical result the published simulations used;
    # the discrepancy is recorded in the vignette Errata.
    # ------------------------------------------------------------------

    e_renalimp_k12_gs704277 <- log(12.33 / 12.53)
    label("Renal impairment on the GS-704277 exchange rate constant (log ratio)")            # Zhang 2025 Table 1 row Q_IM: 12.53 /h healthy, 12.33 /h severe renal impairment

    e_renalimp_k12_gs441524 <- log(0.038 / 0.039)
    label("Renal impairment on the GS-441524 exchange rate constant (log ratio)")            # Zhang 2025 Table 1 row Q_NUC: 0.039 /h healthy, 0.038 /h severe renal impairment

    e_renalimp_kmet_gs704277_central <- log(0.19 / 0.22)
    label("Renal impairment on central remdesivir-to-GS-704277 conversion (log ratio)")      # Zhang 2025 Table 1 row K_C,IM: 0.22 /h healthy, 0.19 /h severe renal impairment. Named in Results 3.1 as covariate-carrying.

    e_renalimp_kmet_gs704277_peripheral1 <- log(0.19 / 0.31)
    label("Renal impairment on peripheral remdesivir-to-GS-704277 conversion (log ratio)")   # Zhang 2025 Table 1 row K_P,IM: 0.31 /h healthy, 0.19 /h severe renal impairment. Named in Results 3.1 as covariate-carrying.

    e_renalimp_kmet_gs441524_central <- log(0.85 / 0.38)
    label("Renal impairment on central GS-704277-to-GS-441524 conversion (log ratio)")       # Zhang 2025 Table 1 row K_C,NUC: 0.38 /h healthy, 0.85 /h severe renal impairment

    e_renalimp_kmet_gs441524_peripheral1 <- log(0.5 / 2.44)
    label("Renal impairment on peripheral GS-704277-to-GS-441524 conversion (log ratio)")    # Zhang 2025 Table 1 row K_P,NUC: 2.44 /h healthy, 0.5 /h severe renal impairment. Named in Results 3.1 as covariate-carrying.

    e_renalimp_kmet_gs443902_peripheral1 <- log(8.08e10 / 3.28e10)
    label("Renal impairment on peripheral GS-441524-to-GS-443902 conversion (log ratio)")    # Zhang 2025 Table 1 row K_P,NTP: 3.28E10 /h healthy, 8.08E10 /h severe renal impairment. Both values act as an instantaneous sink, so the ratio has no discernible effect on the simulated profiles.

    # ------------------------------------------------------------------
    # NO INTER-INDIVIDUAL VARIABILITY IS CARRIED.
    #
    # Zhang 2025 Methods 2.3 states "Random effects were assigned to all
    # parameters with lognormal distributions", but no omega, variance, SD
    # or CV% is reported for any parameter anywhere in the paper, and there
    # is no supplement (the EuropePMC supplementaryFiles bundle for
    # PMC11982744 contains only the six publisher figure files). Inventing
    # variances is not permitted, so this model is typical-value only.
    # The etas are OMITTED rather than written as `~ fixed(0)` because a
    # zero-variance diagonal makes OMEGA singular and breaks the Cholesky
    # sampler used by rxSolve (same handling as
    # Kim_2026_midazolam_postecmo and Thoueille_2026_salmeterol).
    # Recorded in the vignette Errata.
    # ------------------------------------------------------------------

    # ------------------------------------------------------------------
    # RESIDUAL ERROR
    #
    # Zhang 2025 Equation 7 declares a combined additive-plus-proportional
    # error model, C_ij = Y_ij + sqrt(a^2 + (b * Y_ij)^2) * e, "where a and
    # b are additive and proportional errors respectively". Neither a nor b
    # is reported for any of the three analytes. The structure is carried
    # here with both terms fixed at zero so that the declared error model
    # is visible and a re-fit can estimate it; simulations from this file
    # are noise-free. Recorded in the vignette Errata.
    # ------------------------------------------------------------------

    addSd <- fixed(0)
    label("Additive residual SD on remdesivir (ng/mL); magnitude not reported")              # Zhang 2025 Equation 7 declares the combined error model but reports no value for a
    propSd <- fixed(0)
    label("Proportional residual SD on remdesivir (fraction); magnitude not reported")       # Zhang 2025 Equation 7 declares the combined error model but reports no value for b

    addSd_gs704277 <- fixed(0)
    label("Additive residual SD on GS-704277 (ng/mL); magnitude not reported")               # Zhang 2025 Equation 7 declares the combined error model but reports no value for a
    propSd_gs704277 <- fixed(0)
    label("Proportional residual SD on GS-704277 (fraction); magnitude not reported")        # Zhang 2025 Equation 7 declares the combined error model but reports no value for b

    addSd_gs441524 <- fixed(0)
    label("Additive residual SD on GS-441524 (ng/mL); magnitude not reported")               # Zhang 2025 Equation 7 declares the combined error model but reports no value for a
    propSd_gs441524 <- fixed(0)
    label("Proportional residual SD on GS-441524 (fraction); magnitude not reported")        # Zhang 2025 Equation 7 declares the combined error model but reports no value for b
  })

  model({
    # ------------------------------------------------------------------
    # Individual parameters. Renal impairment enters as a multiplicative
    # shift on the log scale, i.e. exp(log(theta_control) + log(ratio) *
    # RENALIMP), so RENALIMP = 0 returns the healthy-control column of
    # Zhang 2025 Table 1 and RENALIMP = 1 returns the severe-renal-
    # impairment column, both exactly as printed.
    # ------------------------------------------------------------------
    k12 <- exp(lk12)
    kel <- exp(lkel)

    k12_gs704277 <- exp(lk12_gs704277 + e_renalimp_k12_gs704277 * RENALIMP)
    k12_gs441524 <- exp(lk12_gs441524 + e_renalimp_k12_gs441524 * RENALIMP)

    kmet_gs704277_central <-
      exp(lkmet_gs704277_central + e_renalimp_kmet_gs704277_central * RENALIMP)
    kmet_gs704277_peripheral1 <-
      exp(lkmet_gs704277_peripheral1 + e_renalimp_kmet_gs704277_peripheral1 * RENALIMP)
    kmet_gs441524_central <-
      exp(lkmet_gs441524_central + e_renalimp_kmet_gs441524_central * RENALIMP)
    kmet_gs441524_peripheral1 <-
      exp(lkmet_gs441524_peripheral1 + e_renalimp_kmet_gs441524_peripheral1 * RENALIMP)
    kmet_gs443902_peripheral1 <-
      exp(lkmet_gs443902_peripheral1 + e_renalimp_kmet_gs443902_peripheral1 * RENALIMP)

    # ------------------------------------------------------------------
    # ODE system -- Zhang 2025 Equations 1-6, transcribed one for one.
    # The states are concentrations, so the exchange terms are written
    # exactly as the paper writes them, Q_i * (C_other - C_this), with no
    # volume ratio. Equation numbers are given against each line.
    #
    #   (1) dC_C,RDV/dt = Q_RDV (C_P,RDV - C_C,RDV) - K_C,IM  C_C,RDV - CL_C,RDV C_C,RDV
    #   (2) dC_P,RDV/dt = Q_RDV (C_C,RDV - C_P,RDV) - K_P,IM  C_P,RDV
    #   (3) dC_C,IM /dt = Q_IM  (C_P,IM  - C_C,IM ) + K_C,IM  C_C,RDV - K_C,NUC C_C,IM
    #   (4) dC_P,IM /dt = Q_IM  (C_C,IM  - C_P,IM ) + K_P,IM  C_P,RDV - K_P,NUC C_P,IM
    #   (5) dC_C,NUC/dt = Q_NUC (C_P,NUC - C_C,NUC) + K_C,NUC C_C,IM
    #   (6) dC_P,NUC/dt = Q_NUC (C_C,NUC - C_P,NUC) + K_P,NUC C_P,IM - K_P,NTP C_P,NUC
    #
    # Note the asymmetry the paper intends and that Figure 2 confirms:
    # remdesivir is eliminated ONLY from its central compartment (no
    # elimination term in equation 2), and GS-441524 is removed ONLY from
    # its peripheral compartment, by conversion to the active triphosphate
    # GS-443902 (no elimination term in equation 5). There is no renal or
    # other clearance term on either metabolite.
    # ------------------------------------------------------------------
    d/dt(central) <-
      k12 * (peripheral1 - central) -
      kmet_gs704277_central * central -
      kel * central                                                          # Equation 1
    d/dt(peripheral1) <-
      k12 * (central - peripheral1) -
      kmet_gs704277_peripheral1 * peripheral1                                # Equation 2

    d/dt(central_gs704277) <-
      k12_gs704277 * (peripheral1_gs704277 - central_gs704277) +
      kmet_gs704277_central * central -
      kmet_gs441524_central * central_gs704277                               # Equation 3
    d/dt(peripheral1_gs704277) <-
      k12_gs704277 * (central_gs704277 - peripheral1_gs704277) +
      kmet_gs704277_peripheral1 * peripheral1 -
      kmet_gs441524_peripheral1 * peripheral1_gs704277                       # Equation 4

    d/dt(central_gs441524) <-
      k12_gs441524 * (peripheral1_gs441524 - central_gs441524) +
      kmet_gs441524_central * central_gs704277                               # Equation 5
    d/dt(peripheral1_gs441524) <-
      k12_gs441524 * (central_gs441524 - peripheral1_gs441524) +
      kmet_gs441524_peripheral1 * peripheral1_gs704277 -
      kmet_gs443902_peripheral1 * peripheral1_gs441524                       # Equation 6

    # ------------------------------------------------------------------
    # Observations. The states already hold plasma concentrations in
    # ng/mL, so no division by a volume is required (and none is possible
    # -- the paper reports no volume). Zhang 2025 Figure 3 plots exactly
    # these three central-compartment quantities against the digitised
    # observations, with ng/mL on every y axis.
    # ------------------------------------------------------------------
    Cc <- central
    Cc_gs704277 <- central_gs704277
    Cc_gs441524 <- central_gs441524

    Cc ~ add(addSd) + prop(propSd)
    Cc_gs704277 ~ add(addSd_gs704277) + prop(propSd_gs704277)
    Cc_gs441524 ~ add(addSd_gs441524) + prop(propSd_gs441524)
  })
}
