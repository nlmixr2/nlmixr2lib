Stevens_2012_remoxipride <- function() {
  description <- "Mechanism-based PK/PD model for the prolactin response to remoxipride in rats: 3-compartment plasma + brain-ECF + peripheral PK with parallel intranasal absorption (systemic and direct nose-to-brain), feeding a pool model for prolactin synthesis (with positive feedback), storage in lactotrophs, release into plasma, and elimination, where remoxipride brain-ECF concentration drives an Emax stimulation of prolactin release"
  reference <- paste(
    "Stevens J., Ploeger B. A., Hammarlund-Udenaes M., Osswald G.,",
    "van der Graaf P. H., Danhof M., de Lange E. C. M. (2012).",
    "Mechanism-based PK-PD model for the prolactin biological system",
    "response following an acute dopamine inhibition challenge:",
    "quantitative extrapolation to humans.",
    "J Pharmacokinet Pharmacodyn 39(5):463-477.",
    "doi:10.1007/s10928-012-9262-4.",
    "DDMORE Foundation Model Repository: DDMODEL00000268.",
    sep = " "
  )
  vignette <- "Stevens_2012_remoxipride"
  units <- list(time = "h", dosing = "mg", concentration = "ug/L (prolactin)")

  ddmore_id <- "DDMODEL00000268"
  replicate_of <- NULL

  covariateData <- list(
    STUDY_DD = list(
      description = "Double-dosing-study indicator: 1 = subject was enrolled in the double-dosing protocol (Stevens 2012 'study 3'), 0 = single-dose (intravenous or intranasal) protocols pooled in the same fit.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (single-dose protocols)",
      notes = "Derived from the source dataset's STUD column (paper-specific integer). The .mod (line 51) sets STDY = 1 when STUD == 3 and 0 otherwise; that mapping is reproduced in covariate notes here. STUDY_DD enters as a multiplicative shift on the typical baseline plasma prolactin (TBSL = bsl * (1 + e_studydd_bsl * STUDY_DD)), making the double-dosing-study baseline ~29% lower than the single-dose-study baseline.",
      source_name = "STUD"
    )
  )

  population <- list(
    n_subjects = NA_integer_,
    n_studies = 3,
    species = "Rattus norvegicus (Wistar rat)",
    weight_range = "approximately 0.20-0.30 kg (per simulated dataset; baseline-prolactin and dose-response cohorts pooled)",
    sex_female_pct = NA_real_,
    disease_state = "healthy male Wistar rats; chronic intracerebral microdialysis cannulae for brain extracellular-fluid (ECF) sampling and femoral-vein cannulae for systemic dosing and plasma sampling",
    dose_range = "Single intravenous remoxipride 4, 8, or 16 mg/kg (dose-response study); double intravenous remoxipride 3.8 mg/kg pulses with varying inter-dose intervals (double-dosing study); a separate baseline-variation study without drug administration. External validation cohort (not contained in the DDMORE bundle's simulated dataset) received 4 / 8 / 16 mg/kg intranasally.",
    notes = "Final estimates from `Output_real_PK_rats.lst` (NONMEM 7.3 FOCEI INTERACTION, OFV 2021.454, MINIMIZATION SUCCESSFUL with a problem-occurred warning, ETAshrinkage 5.3% on the BSL eta with the STDM and K67 etas FIX 0). The PK structural model is fixed from a previously published rat remoxipride PK study (Westerhout et al. 2011, doi:10.1007/s11095-011-0395-8) and reproduced here verbatim from THETA(3..13) FIX in the DDMORE control stream. The PD parameters (TH1, TH2, TH14-TH19) are estimated under an NWPRI prior on TH1/TH2. The model is preclinical (rat-only); no human PK or human PD parameters are encoded in this implementation."
  )

  ini({
    # All values are FINAL ESTIMATES from
    # Output_real_PK_rats.lst (FINAL PARAMETER ESTIMATE block, lines 1283-1284)
    # mapped against the THETA / OMEGA / SIGMA labels in
    # Executable_PK_rats.txt $PK / $THETA / $OMEGA / $SIGMA.
    #
    # PK structural parameters (THETA(3)..THETA(13), all FIX in the .mod
    # because they were carried in from the upstream rat PK paper
    # Westerhout 2011). Final-estimate values match the FIX initial values.

    lcl       <- fixed(log(1.12));      label("Plasma clearance CL3 (L/h)")                                 # TH3 FIX, .lst L1283
    lvc       <- fixed(log(0.0881));    label("Central (plasma) volume V3 (L/kg)")                          # TH4 FIX
    lq_brain  <- fixed(log(0.700));     label("Plasma <-> brain-ECF intercompartmental clearance Q4 (L/h)") # TH5 FIX
    lvbrain   <- fixed(log(0.873));     label("Brain-ECF (microdialysis) volume V4 (L/kg)")                 # TH6 FIX
    lq        <- fixed(log(1.20));      label("Plasma <-> peripheral intercompartmental clearance Q5 (L/h)") # TH7 FIX
    lvp       <- fixed(log(0.417));     label("Peripheral volume V5 (L/kg)")                                # TH8 FIX
    lfk40     <- fixed(log(0.302));     label("FK40 = K40/K30 ratio (brain-ECF outflow rate as a fraction of plasma elimination rate, unitless)") # TH9 FIX
    lka       <- fixed(log(1.54));      label("Intranasal absorption rate KA into plasma (1/h)")            # TH10 FIX
    lftot     <- fixed(log(0.892));     label("Total intranasal bioavailability FTOT (unitless)")           # TH11 FIX
    f_systemic <- fixed(0.249);         label("Fraction of FTOT delivered via the systemic absorption depot (depot -> central); the remaining 1 - f_systemic is delivered via the direct nose-to-brain depot (depot_brain -> brain_csf)") # TH12 FIX
    lk24      <- fixed(log(33.1 / 1000)); label("Direct nose-to-brain absorption rate K24 (1/h); the .mod stores K24 * 1000 in THETA(13) and divides by 1000 inside $PK") # TH13 FIX

    # PD parameters (TH1, TH2, TH14..TH19); estimated under an NWPRI prior on
    # TH1 (kel_prl) and TH2 (ec50_prl). Final-estimate values from .lst.

    lkel_prl    <- log(5.20);   label("Prolactin elimination rate kel_prl = K70 (1/h)")                            # TH1, .lst L1283
    lec50_prl   <- log(0.0510); label("Remoxipride brain-ECF EC50 for stimulation of prolactin release ec50_prl (mg.kg/L; the source's brain-ECF 'concentration' is amount/V_brain with V_brain in L/kg, so the unit carries a kg)") # TH2, .lst L1283
    lbsl        <- log(6.64);   label("Typical baseline plasma prolactin BSL (ug/L)")                              # TH14, .lst L1284
    lemax_prl   <- log(17.5);   label("Maximum stimulation factor of prolactin release Emax_prl (unitless multiplier on the lactotroph release rate)") # TH15, .lst L1284
    lstdm       <- log(0.125);  label("Positive-feedback synthesis-stimulation slope stdm (1/(ug/L))")             # TH16, .lst L1284
    e_studydd_bsl <- -0.290;    label("Double-dosing-study effect on the typical baseline prolactin (additive on (1 + e * STUDY_DD); reference is single-dose)") # TH17, .lst L1284
    lk_release  <- log(0.740);  label("Lactotroph -> plasma prolactin release rate k_release = K67 (1/h)")         # TH18, .lst L1284
    gamma_prl   <- fixed(1.00); label("Hill coefficient on the remoxipride brain-ECF stimulation of prolactin release (FIX 1 -> hyperbolic Emax)") # TH19 FIX, .lst L1284

    # IIV. The .mod fits ETA(1..3) but $OMEGA fixes ETA(2) and ETA(3) to 0
    # in the production run (only the BSL eta is non-zero). Variance from
    # the .lst OMEGA block, line 1294.
    etalbsl ~ 0.0465  # OMEGA(1,1), .lst L1294: log-scale variance on baseline prolactin

    # Residual error. NONMEM .mod $ERROR is
    #   Y = IPRED * (1 + EPS(1)) + EPS(2)
    # so EPS(1) is proportional and EPS(2) is additive (in ug/L). SIGMA
    # variances are 0.0780 (prop) and 6.62 (add) per .lst lines 1310, 1313;
    # nlmixr2 uses SD on the same scale.
    propSd <- sqrt(0.0780); label("Proportional residual error on plasma prolactin (fraction)") # SIGMA(1,1), .lst L1310 -> sqrt(variance)
    addSd  <- sqrt(6.62);   label("Additive residual error on plasma prolactin (ug/L)")         # SIGMA(2,2), .lst L1313 -> sqrt(variance)
  })

  model({
    # =====================================================================
    # Verbatim translation of Executable_PK_rats.txt $PK + $DES + $ERROR.
    # NONMEM A(i) compartments map to nlmixr2 named states as follows:
    #   A(1) = depot          (intranasal "fast" depot, rate KA -> central)
    #   A(2) = depot_brain    (intranasal "slow" depot, rate K24 -> brain_csf)
    #   A(3) = central        (plasma, DEFDOSE for IV)
    #   A(4) = brain_csf      (microdialysis brain extracellular fluid)
    #   A(5) = peripheral1    (PK peripheral)
    #   A(6) = lactotroph     (prolactin storage in lactotrophs)
    #   A(7) = prolactin      (plasma prolactin, the observed)
    # =====================================================================

    # Individual PK parameters (all log-transformed and FIX, no etas).
    cl       <- exp(lcl)
    vc       <- exp(lvc)
    q_brain  <- exp(lq_brain)
    vbrain   <- exp(lvbrain)
    q        <- exp(lq)
    vp       <- exp(lvp)
    fk40     <- exp(lfk40)
    ka       <- exp(lka)
    k24      <- exp(lk24)
    ftot     <- exp(lftot)

    # Individual PD parameters; etalbsl is the only non-FIX random effect.
    bsl_typ  <- exp(lbsl) * (1 + e_studydd_bsl * STUDY_DD)        # .mod L56: TBSL = THETA(14) * (1 + THETA(17) * STDY)
    bsl_i    <- bsl_typ * exp(etalbsl)                            # .mod L57: BSL = TBSL * EXP(ETA(1))
    kel_prl  <- exp(lkel_prl)                                     # .mod L60: K70 = THETA(1)
    ec50_prl <- exp(lec50_prl)                                    # .mod L62
    emax_prl <- exp(lemax_prl)                                    # .mod L61
    stdm     <- exp(lstdm)                                        # .mod L63 (ETA(2) FIX 0)
    k_release <- exp(lk_release)                                  # .mod L59 (ETA(3) FIX 0)

    # Micro-rate constants derived from CL/V parameterization (.mod L40-45).
    kel_plasma <- cl / vc                                         # K30 = CL3 / V3
    k40        <- fk40 * kel_plasma                               # K40 = FK40 * K30
    k34        <- q_brain / vc                                    # K34 = Q4 / V3
    k43        <- q_brain / vbrain                                # K43 = Q4 / V4
    k35        <- q / vc                                          # K35 = Q5 / V3
    k53        <- q / vp                                          # K53 = Q5 / V5

    # Initial conditions. The .mod sets baseline lactotroph and plasma
    # prolactin amounts so that the system is at steady state in the
    # absence of drug:
    #   A_0(6) = (BSL * K70) / K67    (.mod L66)
    #   A_0(7) = BSL                  (.mod L67)
    lactotroph(0) <- bsl_i * kel_prl / k_release
    prolactin(0)  <- bsl_i

    # =====================================================================
    # ODE system (.mod $DES, lines 74-99)
    # =====================================================================
    d/dt(depot)       <- -ka * depot                                                # .mod L75
    d/dt(depot_brain) <- -k24 * depot_brain                                         # .mod L76
    d/dt(central)     <-  ka * depot - kel_plasma * central -
                          k34 * central - k35 * central +
                          k43 * brain_csf + k53 * peripheral1                       # .mod L77
    d/dt(brain_csf)   <-  k24 * depot_brain + k34 * central -
                          k43 * brain_csf - k40 * brain_csf                         # .mod L78
    d/dt(peripheral1) <-  k35 * central - k53 * peripheral1                         # .mod L88

    # Bioavailability split for intranasal dosing. The .mod (L36-38) divides
    # the total intranasal bioavailability FTOT into a systemic-route
    # fraction (F1 = FTOT * THETA(12)) and a direct-nose-to-brain fraction
    # (F2 = FTOT - F1). Intravenous doses bypass both depots and go directly
    # to central; users dose IV by setting the cmt column to `central`.
    f(depot)       <- ftot * f_systemic
    f(depot_brain) <- ftot * (1 - f_systemic)

    # Drug effect on prolactin release. CP is the "concentration" of
    # remoxipride in brain ECF; in the source it is a body-weight-
    # normalized amount-per-(L/kg) since V4 is L/kg, so CP carries an
    # extra factor of kg (mg.kg/L). EC50 is on the same scale.
    cp_brain <- brain_csf / vbrain
    eff_num  <- emax_prl * cp_brain^gamma_prl
    eff_den  <- ec50_prl^gamma_prl + cp_brain^gamma_prl
    eff      <- 1 + eff_num / (eff_den + 1e-30)                                     # .mod L82

    # Positive feedback floor: the synthesis stimulation uses A(7) but
    # clamped from below at the baseline so the term stdm * (lna7 - bsl)
    # is non-negative (.mod L85-93).
    lna7 <- prolactin
    if (prolactin < bsl_i) lna7 <- bsl_i

    std <- stdm * (lna7 - bsl_i)                                                    # .mod L95
    kf  <- bsl_i * kel_prl                                                          # .mod L68: KF = BSL * K70
    kft <- kf * (1 + std)                                                           # .mod L96

    d/dt(lactotroph) <- kft - lactotroph * k_release * eff                          # .mod L97
    d/dt(prolactin)  <- lactotroph * k_release * eff - prolactin * kel_prl          # .mod L98

    # Observation. The .mod $ERROR (L102-103) is
    #   Y = IPRED * (1 + EPS(1)) + EPS(2)
    # with IPRED = A(7), so this is a combined proportional + additive
    # error on plasma prolactin in ug/L.
    prolactin ~ prop(propSd) + add(addSd)
  })
}
