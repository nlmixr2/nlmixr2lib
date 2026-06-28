Westerhout_2013_quinidine <- function() {
  description <- paste(
    "Preclinical (rat, male Wistar WU). Systems-based pharmacokinetic",
    "(SBPK) model for quinidine intra-brain distribution following IV",
    "infusion, fit jointly to unbound plasma, brain parenchymal",
    "extracellular fluid (brain ECF), CSF in the lateral ventricle",
    "(CSF_LV) and cisterna magna (CSF_CM), and end-of-experiment total",
    "(deep) brain concentrations, with simultaneous mechanistic state",
    "compartments for CSF in the combined third + fourth ventricles",
    "(CSF_TFV) and the subarachnoid space (CSF_SAS) carrying the",
    "ventricular CSF flow from LV through SAS back to systemic plasma",
    "at the fixed physiological rate Q_CSF = 2.2 uL/min, plus a",
    "brain-ECF-to-CSF_LV flow at Q_ECF = 0.2 uL/min (Westerhout 2013,",
    "J Pharmacokinet Pharmacodyn). Two systemic peripheral",
    "compartments (V_PER1, V_PER2) with inter-compartmental",
    "clearances Q_PL-PER1, Q_PL-PER2 carry the multi-exponential",
    "plasma decline. P-glycoprotein (P-gp) activity (binary indicator",
    "CONMED_TARIQUIDAR = 0 control / 1 tariquidar-inhibited)",
    "modulates the systemic elimination CL_E (1.9-fold increase when",
    "P-gp is active per Table 4) and every transfer clearance between",
    "plasma and the brain compartments: passive influx into each",
    "brain compartment is reduced when P-gp is active (influx",
    "hindrance, subtractive) and passive efflux is increased (efflux",
    "enhancement, additive). At the BCSFB level, P-gp acts as an",
    "efflux transporter at the LV (CL_LV-PL,P-gp estimated to 0 in",
    "the combined model, with the P-gp-mediated component carried by",
    "CL_PL-LV,P-gp) and is absent at the CM (both P-gp components",
    "fixed to 0 in the combined model). The plasma-to-TFV transfer",
    "clearance is structurally assumed equal to plasma-to-LV (no TFV",
    "microdialysis sampling). Parameter values are the paper's",
    "preferred 'efflux enhancement + influx hindrance' (combined) SBPK",
    "model from Table 4, column 3 (OFV = 17,969); the two alternative",
    "P-gp mechanism variants (efflux-enhancement-only, OFV = 18,105;",
    "influx-hindrance-only, OFV = 18,030) and the simpler preliminary",
    "compartmental model (Table 3) are discussed in the validation",
    "vignette but not extracted as separate model files."
  )
  reference <- paste(
    "Westerhout J, Smeets J, Danhof M, de Lange ECM.",
    "The impact of P-gp functionality on non-steady state",
    "relationships between CSF and brain extracellular fluid.",
    "J Pharmacokinet Pharmacodyn. 2013;40(3):327-342.",
    "doi:10.1007/s10928-013-9314-4."
  )
  vignette <- "Westerhout_2013_quinidine"

  units <- list(
    time          = "min",
    dosing        = "ng",
    concentration = paste(
      "ng/mL (unbound plasma quinidine after correction for plasma",
      "protein binding; unbound brain ECF / CSF concentrations after",
      "correction for in vivo microdialysis probe recovery; deep brain",
      "concentrations are total brain after subtraction of the brain",
      "ECF contribution)"
    )
  )

  covariateData <- list(
    CONMED_TARIQUIDAR = list(
      description        = paste(
        "Indicator for tariquidar (XR9576) pre-administration:",
        "1 = rat received a 15 mg/kg IV tariquidar infusion in 5%",
        "glucose / saline vehicle (100 uL/min/kg over 10 min) 30 min",
        "before the quinidine infusion; 0 = rat received vehicle alone.",
        "Tariquidar is a selective third-generation P-glycoprotein",
        "(P-gp) inhibitor; in this design it is assumed to fully",
        "inhibit P-gp at the BBB, BCSFB, and intra-parenchymal sites",
        "throughout the 6 h experimental window (paper Discussion,",
        "page 339: '50-fold higher than the IC50 value up to 3 h after",
        "administration')."
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (vehicle co-administration; no tariquidar; P-gp fully active)",
      notes              = paste(
        "Time-fixed per animal in Westerhout 2013. Source paper",
        "Materials and Methods 'Experimental set-up' (page 329)",
        "describes the tariquidar dosing protocol; Methods 'Animals'",
        "(page 328) details the 4-arm 2 x 2 factorial design",
        "(quinidine dose 10 or 20 mg/kg crossed with TQD- or TQD+).",
        "The covariate enters the model multiplicatively through the",
        "P-gp-mediated component of every transfer clearance: when",
        "CONMED_TARIQUIDAR = 1 the P-gp components are switched off",
        "and the effective clearances collapse to their passive",
        "values; when CONMED_TARIQUIDAR = 0 the P-gp components add",
        "(efflux) or subtract (influx hindrance) per the Westerhout",
        "Appendix mass-balance equations."
      ),
      source_name        = "TARIQUIDAR"
    )
  )

  population <- list(
    species        = "rat (male Wistar WU, Charles River Maastricht NL)",
    n_subjects     = 48L,
    n_studies      = 1L,
    age_range      = "adult (not reported in days; rats arrived weighing 225-275 g and were housed >= 5-7 days then individually housed for 7 days post-surgery)",
    weight_range   = "225-275 g body weight on arrival",
    sex_female_pct = 0,
    race_ethnicity = NA,
    disease_state  = paste(
      "Healthy adult male Wistar WU rats with chronically implanted",
      "femoral arterial and venous cannulae (blood sampling and IV",
      "drug administration) and two CMA/12 microdialysis guides in",
      "combinations of striatum (ST, for brain ECF), lateral",
      "ventricle (LV, for CSF_LV), and cisterna magna (CM, for",
      "CSF_CM). End-of-experiment whole-brain isolation supplied",
      "total brain concentrations used to back-calculate deep brain",
      "(brain intracellular) concentrations."
    ),
    dose_range     = paste(
      "Single IV quinidine infusion of 10 or 20 mg/kg in saline",
      "delivered at 100 uL/min/kg over 10 min (e.g., 250 ug/min over",
      "10 min in a 250 g rat for the 10 mg/kg arm), with infusion",
      "start corrected for internal tubing volume so that t = 0 min",
      "is the actual infusion start. Tariquidar pre-administration",
      "arm: 15 mg/kg IV in 5% glucose / saline at 100 uL/min/kg over",
      "10 min, 30 min before quinidine."
    ),
    regions        = "preclinical (in-vivo rat); Leiden / Amsterdam Center for Drug Research, The Netherlands",
    notes          = paste(
      "60 male Wistar WU rats total: 12 dedicated to the in-vivo",
      "microdialysis probe recovery determination (reverse",
      "dialysis, retrodialysis calibration), 48 to the brain",
      "disposition fits (n = 12 per cell of the 2 x 2 factorial:",
      "10 mg/kg + vehicle, 10 mg/kg + tariquidar, 20 mg/kg + vehicle,",
      "20 mg/kg + tariquidar). Plasma protein binding of quinidine",
      "was independently measured at 86.5 +/- 5.5% (linear, not",
      "tariquidar-dependent); only unbound plasma concentrations",
      "enter the model (the model state `central` is the unbound",
      "quinidine plasma amount A_pl,u in the paper's notation).",
      "Microdialysate concentrations from ST, LV, and CM probes were",
      "corrected for the location-specific in vivo recoveries (9.1%,",
      "2.9%, 3.5% respectively, all unchanged by tariquidar) before",
      "entering the model."
    )
  )

  ini({
    # ------------------------------------------------------------
    # Systemic elimination. Paper Table 4, SBPK 'efflux enhancement
    # and influx hindrance' column (the paper's preferred final
    # model, OFV = 17,969). Appendix Eq.: k_E = (CL_E,p + CL_E,P-gp)
    # / V_PL. CL_E,p = 95.9 mL/min is reported directly; the additive
    # P-gp contribution CL_E,P-gp is not tabulated as a column but is
    # implied by the reported "P-gp effect on CL_E = 1.9-fold
    # increase", giving CL_E,P-gp = (1.9 - 1) * 95.9 = 86.31 mL/min.
    # The additive split is the paper's structural form; using a
    # multiplicative covariate effect on a single CL_E baseline would
    # also reproduce the typical-value behaviour but would not match
    # the Westerhout Appendix equations.
    lcl     <- log(95.9);   label("Passive systemic elimination CL_E,p (mL/min) - P-gp inhibited (tariquidar present)")  # Table 4, SBPK combined: 95.9 +/- 11.0 mL/min
    lcl_pgp <- log(86.31);  label("Additive P-gp-mediated elimination CL_E,P-gp (mL/min) - active only when CONMED_TARIQUIDAR = 0")  # Derived from Table 4 SBPK combined: P-gp effect on CL_E = 1.9 +/- 0.2-fold increase; (1.9 - 1) * 95.9 = 86.31 mL/min

    # ------------------------------------------------------------
    # Inter-compartmental plasma <-> peripheral.
    lq      <- log(1190);   label("Plasma <-> peripheral1 inter-compartmental clearance Q_PL-PER1 (mL/min)")  # Table 4: 1190 +/- 135 mL/min
    lq2     <- log(333);    label("Plasma <-> peripheral2 inter-compartmental clearance Q_PL-PER2 (mL/min)")  # Table 4: 333 +/- 94 mL/min

    # ------------------------------------------------------------
    # Volumes (mL). V_PL is fixed at the Davies 1993 literature
    # value (paper Table 4 cites it via reference [52]); V_PER1 and
    # V_PER2 are estimated.
    lvc     <- fixed(log(10.6));   label("Plasma volume V_PL (mL), fixed at literature value")  # Table 4: 10.6 mL [52]
    lvp     <- log(6800);          label("First peripheral volume V_PER1 (mL)")                  # Table 4: 6.8 +/- 1.7 L = 6800 mL
    lvp2    <- log(13300);         label("Second peripheral volume V_PER2 (mL)")                 # Table 4: 13.3 +/- 2.2 L = 13300 mL

    # ------------------------------------------------------------
    # Plasma <-> brain_deep transport.
    # Paper Appendix:
    #   k_PL-DBR = (CL_PL-DBR,p - CL_PL-DBR,P-gp) / V_PL   (influx; P-gp HINDERS)
    #   k_DBR-PL = (CL_DBR-PL,p + CL_DBR-PL,P-gp) / V_DBR  (efflux; P-gp ENHANCES)
    # Table 4 values in uL/min; encoded here in mL/min (uL/min / 1000).
    lcl_pl_dbr_p   <- log(2.180);   label("Passive plasma -> brain_deep clearance CL_PL-DBR,p (mL/min) - tariquidar value")           # Table 4: 2180 +/- 384 uL/min
    lcl_pl_dbr_pgp <- log(1.900);   label("P-gp influx hindrance plasma -> brain_deep CL_PL-DBR,P-gp (mL/min) - subtracted when P-gp active")  # Table 4: 1900 +/- 373 uL/min
    lcl_dbr_pl_p   <- log(0.0372);  label("Passive brain_deep -> plasma clearance CL_DBR-PL,p (mL/min)")                              # Table 4: 37.2 +/- 7.2 uL/min
    lcl_dbr_pl_pgp <- log(0.0196);  label("P-gp efflux enhancement brain_deep -> plasma CL_DBR-PL,P-gp (mL/min) - added when P-gp active")  # Table 4: 19.6 +/- 10.9 uL/min

    # ------------------------------------------------------------
    # Plasma <-> brain_ecf transport.
    lcl_pl_ecf_p   <- log(0.0502);  label("Passive plasma -> brain_ecf clearance CL_PL-ECF,p (mL/min)")                               # Table 4: 50.2 +/- 5.0 uL/min
    lcl_pl_ecf_pgp <- log(0.0338);  label("P-gp influx hindrance plasma -> brain_ecf CL_PL-ECF,P-gp (mL/min)")                        # Table 4: 33.8 +/- 5.1 uL/min
    lcl_ecf_pl_p   <- log(0.0063);  label("Passive brain_ecf -> plasma clearance CL_ECF-PL,p (mL/min)")                               # Table 4: 6.3 +/- 0.8 uL/min
    lcl_ecf_pl_pgp <- log(0.0053);  label("P-gp efflux enhancement brain_ecf -> plasma CL_ECF-PL,P-gp (mL/min)")                      # Table 4: 5.3 +/- 1.7 uL/min

    # ------------------------------------------------------------
    # Plasma <-> brain_csf_lv transport (and identically plasma <->
    # brain_csf_tfv per the paper's structural assumption that the
    # plasma-to-TFV transfer clearance equals the plasma-to-LV
    # transfer clearance, because no TFV microdialysis was performed
    # to identify it separately - paper Methods 'Systems-based
    # modeling approach', page 335).
    lcl_pl_lv_p    <- log(0.009);    label("Passive plasma -> CSF_LV clearance CL_PL-LV,p (mL/min) - also applied to plasma -> CSF_TFV")  # Table 4: 9.0 +/- 0.9 uL/min
    lcl_pl_lv_pgp  <- log(0.0038);   label("P-gp influx hindrance plasma -> CSF_LV CL_PL-LV,P-gp (mL/min) - also applied to plasma -> CSF_TFV")  # Table 4: 3.8 +/- 0.8 uL/min
    lcl_lv_pl_p    <- log(0.00004);  label("Passive CSF_LV -> plasma clearance CL_LV-PL,p (mL/min) - also applied to CSF_TFV -> plasma")    # Table 4: 0.04 +/- 0.01 uL/min
    # CL_LV-PL,P-gp = 0 in the SBPK combined model (Table 4; paper
    # Results 'Modeling P-gp-mediated transport', page 336: "the data
    # were best described by a model with P-gp as an efflux
    # transporter at the BCSFB for LV", encoded via the CL_PL-LV,P-gp
    # component above with the LV-to-plasma P-gp component pinned to
    # 0). Omitted from ini() because log(0) is undefined; the
    # corresponding effective efflux in model() drops the P-gp term.

    # ------------------------------------------------------------
    # Plasma <-> brain_csf_cm transport. In the SBPK combined model
    # (Table 4) both P-gp components for the CM are fixed to 0 (paper
    # Results, page 336: "P-gp ... is absent at the CM"). Only the
    # passive components are estimated.
    lcl_pl_cm_p    <- log(0.0011);   label("Passive plasma -> CSF_CM clearance CL_PL-CM,p (mL/min) - P-gp components fixed to 0")           # Table 4: 1.1 +/- 0.3 uL/min
    lcl_cm_pl_p    <- log(0.0041);   label("Passive CSF_CM -> plasma clearance CL_CM-PL,p (mL/min)")                                        # Table 4: 4.1 +/- 0.5 uL/min

    # ------------------------------------------------------------
    # Inter-animal variability. Westerhout 2013 reports a single eta
    # in Table 4 (eta_CL10 = 0.14 +/- 0.06 in the SBPK combined
    # column), interpreted as the variance of the log-normal random
    # effect on the effective systemic elimination CL_E. The IIV is
    # multiplicative on the additively-combined (passive + P-gp)
    # effective CL_E, matching the NONMEM convention for an
    # additively-decomposed THETA expression with a single ETA on the
    # combined value.
    etalcl   ~ 0.14   # Table 4 SBPK combined: eta_CL10 = 0.14 +/- 0.06 (variance on log-normal IIV)

    # ------------------------------------------------------------
    # Residual error: independent proportional models per observation
    # stream. Westerhout 2013 Table 4 reports the per-output residual
    # variance epsilon (NONMEM $SIGMA convention: y_obs = y_pred *
    # (1 + eps), eps ~ N(0, sigma^2)); propSd is entered as
    # sqrt(sigma^2) so that ~prop(propSd) interprets it as the SD on
    # the multiplicative noise term.
    propSd                <- sqrt(0.22);  label("Proportional residual SD on unbound plasma Cc (fraction)")                                    # Table 4 SBPK combined: eps_PL  = 0.22 +/- 0.03 -> SD 0.469
    propSd_Cbrain_deep    <- sqrt(0.07);  label("Proportional residual SD on deep brain Cbrain_deep (fraction)")                               # Table 4 SBPK combined: eps_DBR = 0.07 +/- 0.02 -> SD 0.265
    propSd_Cbrain_ecf     <- sqrt(0.06);  label("Proportional residual SD on brain ECF Cbrain_ecf (fraction)")                                 # Table 4 SBPK combined: eps_ECF = 0.06 +/- 0.01 -> SD 0.245
    propSd_Cbrain_csf_lv  <- sqrt(0.11);  label("Proportional residual SD on CSF lateral-ventricle Cbrain_csf_lv (fraction)")                  # Table 4 SBPK combined: eps_LV  = 0.11 +/- 0.02 -> SD 0.332
    propSd_Cbrain_csf_cm  <- sqrt(0.08);  label("Proportional residual SD on CSF cisterna-magna Cbrain_csf_cm (fraction)")                     # Table 4 SBPK combined: eps_CM  = 0.08 +/- 0.02 -> SD 0.283
  })

  model({
    # ----- Tariquidar / P-gp indicator -----
    # CONMED_TARIQUIDAR = 1 when tariquidar is co-administered
    # (P-gp inhibited); pgp_active is the complement so that the
    # P-gp-mediated clearance components contribute only in control
    # animals (CONMED_TARIQUIDAR = 0).
    pgp_active <- 1 - CONMED_TARIQUIDAR

    # ----- Volumes (mL) -----
    vc   <- exp(lvc)
    vp   <- exp(lvp)
    vp2  <- exp(lvp2)

    # ----- Brain compartment volumes (mL; literature-fixed rat
    # physiology, Westerhout 2013 Methods 'PK data analysis' page 329
    # and Table 4 footnotes) -----
    v_dbr <- 1.44      # mL - deep (intracellular) brain volume V_DBR (ref [38])
    v_ecf <- 0.290     # mL - brain extracellular fluid volume V_ECF, 290 uL (ref [39])
    v_lv  <- 0.050     # mL - CSF lateral-ventricle volume V_LV, 50 uL (refs [41, 42])
    v_tfv <- 0.050     # mL - CSF third + fourth ventricle volume V_TFV, 50 uL (ref [43])
    v_cm  <- 0.017     # mL - CSF cisterna-magna volume V_CM, 17 uL (refs [44, 45])
    v_sas <- 0.180     # mL - CSF subarachnoid-space volume V_SAS, 180 uL (refs [40, 43])

    # ----- Physiological flow rates (mL/min; literature-fixed
    # values from Westerhout 2013 Table 4) -----
    q_ecf <- 0.0002    # mL/min - brain ECF flow Q_ECF, 0.2 uL/min (refs [39, 46])
    q_csf <- 0.0022    # mL/min - CSF flow Q_CSF, 2.2 uL/min (ref [47])

    # ----- Effective elimination CL_E with additive P-gp component
    # (paper Appendix: k_E = (CL_E,p + CL_E,P-gp) / V_PL) and a single
    # multiplicative log-normal IIV on the effective CL_E -----
    cl <- (exp(lcl) + pgp_active * exp(lcl_pgp)) * exp(etalcl)

    # ----- Effective transfer clearances (mL/min) with P-gp split.
    # Influx: P-gp HINDERS, so effective = passive - P-gp.
    # Efflux: P-gp ENHANCES, so effective = passive + P-gp.
    # Author structural assumption: plasma <-> CSF_TFV uses the same
    # passive and P-gp clearances as plasma <-> CSF_LV (no TFV
    # microdialysis identifies separate values).
    cl_pl_dbr <- exp(lcl_pl_dbr_p) - pgp_active * exp(lcl_pl_dbr_pgp)
    cl_dbr_pl <- exp(lcl_dbr_pl_p) + pgp_active * exp(lcl_dbr_pl_pgp)
    cl_pl_ecf <- exp(lcl_pl_ecf_p) - pgp_active * exp(lcl_pl_ecf_pgp)
    cl_ecf_pl <- exp(lcl_ecf_pl_p) + pgp_active * exp(lcl_ecf_pl_pgp)
    cl_pl_lv  <- exp(lcl_pl_lv_p)  - pgp_active * exp(lcl_pl_lv_pgp)
    cl_lv_pl  <- exp(lcl_lv_pl_p)                                       # CL_LV-PL,P-gp = 0 in SBPK combined
    cl_pl_tfv <- cl_pl_lv                                               # structural: plasma <-> CSF_TFV = plasma <-> CSF_LV
    cl_tfv_pl <- cl_lv_pl                                               # structural: CSF_TFV <-> plasma = CSF_LV <-> plasma
    cl_pl_cm  <- exp(lcl_pl_cm_p)                                       # CL_PL-CM,P-gp = 0 in SBPK combined
    cl_cm_pl  <- exp(lcl_cm_pl_p)                                       # CL_CM-PL,P-gp = 0 in SBPK combined

    # ----- Mass-balance rate constants (1/min). Each k_<from>_<to>
    # is the corresponding effective clearance divided by the SOURCE
    # compartment volume, per the paper Appendix definitions. -----
    kE         <- cl       / vc
    k_pl_per1  <- exp(lq)  / vc
    k_per1_pl  <- exp(lq)  / vp
    k_pl_per2  <- exp(lq2) / vc
    k_per2_pl  <- exp(lq2) / vp2
    k_pl_dbr   <- cl_pl_dbr / vc
    k_dbr_pl   <- cl_dbr_pl / v_dbr
    k_pl_ecf   <- cl_pl_ecf / vc
    k_ecf_pl   <- cl_ecf_pl / v_ecf
    k_pl_lv    <- cl_pl_lv  / vc
    k_lv_pl    <- cl_lv_pl  / v_lv
    k_pl_tfv   <- cl_pl_tfv / vc
    k_tfv_pl   <- cl_tfv_pl / v_tfv
    k_pl_cm    <- cl_pl_cm  / vc
    k_cm_pl    <- cl_cm_pl  / v_cm

    # ----- ODE system (Westerhout 2013 Appendix mass-balance form).
    # The state `central` is the unbound plasma amount A_pl,u; the IV
    # quinidine infusion is dosed via the event table on the central
    # compartment after conversion of the total dose to the unbound
    # fraction (fu_plasma = 0.135 = 1 - 0.865, paper Results page 332).
    # The brain compartments hold drug amount in ng; concentrations
    # are derived by dividing by the corresponding volume. -----
    d/dt(central)       <- -k_pl_per1 * central + k_per1_pl * peripheral1 -
                            k_pl_per2 * central + k_per2_pl * peripheral2 -
                            k_pl_dbr  * central + k_dbr_pl  * brain_deep -
                            k_pl_ecf  * central + k_ecf_pl  * brain_ecf -
                            k_pl_lv   * central + k_lv_pl   * brain_csf_lv -
                            k_pl_tfv  * central + k_tfv_pl  * brain_csf_tfv -
                            k_pl_cm   * central + k_cm_pl   * brain_csf_cm +
                            (q_csf / v_sas) * brain_csf_sas - kE * central
    d/dt(peripheral1)   <- k_pl_per1 * central - k_per1_pl * peripheral1
    d/dt(peripheral2)   <- k_pl_per2 * central - k_per2_pl * peripheral2
    d/dt(brain_deep)    <- k_pl_dbr  * central - k_dbr_pl  * brain_deep
    d/dt(brain_ecf)     <- k_pl_ecf  * central - k_ecf_pl  * brain_ecf -
                            (q_ecf / v_ecf) * brain_ecf
    d/dt(brain_csf_lv)  <- k_pl_lv   * central - k_lv_pl   * brain_csf_lv +
                            (q_ecf / v_ecf) * brain_ecf -
                            (q_csf / v_lv)  * brain_csf_lv
    d/dt(brain_csf_tfv) <- k_pl_tfv  * central - k_tfv_pl  * brain_csf_tfv +
                            (q_csf / v_lv)  * brain_csf_lv -
                            (q_csf / v_tfv) * brain_csf_tfv
    d/dt(brain_csf_cm)  <- k_pl_cm   * central - k_cm_pl   * brain_csf_cm +
                            (q_csf / v_tfv) * brain_csf_tfv -
                            (q_csf / v_cm)  * brain_csf_cm
    d/dt(brain_csf_sas) <- (q_csf / v_cm)  * brain_csf_cm -
                            (q_csf / v_sas) * brain_csf_sas

    # ----- Algebraic concentrations -----
    Cc             <- central       / vc       # unbound plasma quinidine (ng/mL)
    Cbrain_deep    <- brain_deep    / v_dbr
    Cbrain_ecf     <- brain_ecf     / v_ecf
    Cbrain_csf_lv  <- brain_csf_lv  / v_lv
    Cbrain_csf_cm  <- brain_csf_cm  / v_cm
    # Cbrain_csf_tfv and Cbrain_csf_sas are mechanistic state
    # concentrations exposed for diagnostic / simulation use; the
    # original Westerhout 2013 dataset has no microdialysis
    # observations from these compartments and therefore no residual
    # error stream is assigned.
    Cbrain_csf_tfv <- brain_csf_tfv / v_tfv
    Cbrain_csf_sas <- brain_csf_sas / v_sas

    # ----- Observations (proportional error per output stream) -----
    Cc             ~ prop(propSd)
    Cbrain_deep    ~ prop(propSd_Cbrain_deep)
    Cbrain_ecf     ~ prop(propSd_Cbrain_ecf)
    Cbrain_csf_lv  ~ prop(propSd_Cbrain_csf_lv)
    Cbrain_csf_cm  ~ prop(propSd_Cbrain_csf_cm)
  })
}
