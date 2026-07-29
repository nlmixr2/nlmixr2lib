Syvanen_2012_quinidine_rat <- function() {
  description <- paste(
    "Preclinical (rat, male Sprague-Dawley).",
    "Two-compartment plasma + brain extracellular fluid (ECF)",
    "population PK model for quinidine in rats with hippocampal",
    "microdialysis sampling, fit jointly to plasma, brain ECF, and",
    "end-of-experiment total brain concentrations (Syvanen 2012).",
    "The plasma 2-cmt system (central V1 / peripheral V2) couples to",
    "a brain ECF compartment V_Br via asymmetric BBB clearances",
    "(Q_in = f1 * Q_out into brain, Q_out out of brain, Q_out FIXED",
    "at 10.8 mL/min from the paper's bootstrap-stability analysis).",
    "A third observed output, total brain tissue concentration",
    "Cbrain_deep, is modelled as an algebraic equilibrium multiple",
    "of brain ECF (Cbrain_deep = f2 * Cbrain_csf) because the paper",
    "could not estimate separate rate constants for the deep brain",
    "compartment. Two binary covariates: CONMED_TARIQUIDAR (15 mg/kg",
    "IP tariquidar pre-administered 30 min before quinidine, a",
    "selective P-glycoprotein inhibitor) modifies CL, Q_out, Q_in,",
    "and f2; DIS_POSTSE_KAINATE (1 week post-kainate-induced status",
    "epilepticus, rat temporal-lobe-epilepsy paradigm) modifies CL,",
    "V2, and V_Br."
  )
  reference <- paste(
    "Syvanen S, Schenke M, van den Berg D-J, Voskuyl RA,",
    "de Lange ECM. Alteration in P-glycoprotein Functionality",
    "Affects Intrabrain Distribution of Quinidine More Than Brain",
    "Entry - A Study in Rats Subjected to Status Epilepticus by",
    "Kainate. AAPS J. 2012;14(1):87-96.",
    "doi:10.1208/s12248-011-9318-1."
  )
  vignette <- "Syvanen_2012_quinidine_rat"
  paper_specific_etas <- c("etalqin", "etalvbr")

  units <- list(
    time          = "min",
    dosing        = "ng",
    concentration = "ng/mL"
  )

  covariateData <- list(
    DIS_POSTSE_KAINATE = list(
      description        = paste(
        "Indicator for chronic post-status-epilepticus (post-SE) state",
        "in the kainate-induced rat temporal-lobe-epilepsy paradigm:",
        "1 = rat previously subjected to chemoconvulsant SE by",
        "repetitive IP kainic acid (initial 10 mg/kg, then 5 mg/kg",
        "every 30-60 min until stage IV/V Racine seizures or a",
        "30 mg/kg cumulative cap) 6-7 days before microdialysis",
        "surgery, sampled at 7 days post-SE during the pre-spontaneous-",
        "epilepsy window in which P-gp expression has been reported to",
        "peak; 0 = naive control or saline-injected sham control rat."
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (naive control or saline-sham control rat)",
      notes              = paste(
        "Time-fixed per animal in Syvanen 2012 (each rat is allocated",
        "to exactly one of {control, kainate-treated}). Source paper",
        "Materials and Methods 'Animals' (page 88) describes the",
        "kainate dosing protocol; Table I lists per-arm rat counts.",
        "The post-SE observation lag is 7 days from SE induction to",
        "quinidine sampling. Modifies CL (0.537-fold), V2 / peripheral",
        "plasma (0.678-fold), and V_Br / brain ECF (0.592-fold) per",
        "Table II."
      ),
      source_name        = "KAINATE"
    ),
    CONMED_TARIQUIDAR = list(
      description        = paste(
        "Indicator for tariquidar (XR9576) pre-administration:",
        "1 = rat received 15 mg/kg IP tariquidar in 5% glucose / saline",
        "vehicle 30 min before the quinidine intravenous infusion;",
        "0 = rat received vehicle (5% glucose in saline) alone at the",
        "same time point. Tariquidar is a selective third-generation",
        "P-glycoprotein inhibitor that blocks BBB P-gp efflux."
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (vehicle co-administration; no tariquidar)",
      notes              = paste(
        "Time-fixed per animal in Syvanen 2012. Source paper Materials",
        "and Methods 'Pharmacokinetic Study in Rats' (page 89)",
        "describes the tariquidar protocol; Table I lists per-arm rat",
        "counts (each control dose arm split 1:1 vehicle / TQD;",
        "kainate arm split 8:8 vehicle / TQD). Modifies CL (0.835-fold),",
        "Q_out (0.366-fold), Q_in (2.65-fold), and f2 (5.66-fold) per",
        "Table II. The combined Q_in / Q_out shift raises brain ECF",
        "concentrations 7.2-fold (2.65 / 0.366) and the combined f2",
        "shift further raises total brain concentrations to about",
        "40-fold higher than baseline (paper Results 'TQD Treatment',",
        "page 91)."
      ),
      source_name        = "TARIQUIDAR"
    )
  )

  population <- list(
    species        = "rat (male Sprague-Dawley, Harlan)",
    n_subjects     = 64L,
    n_studies      = 1L,
    age_range      = "adult (not reported in days; rats arrived weighing 200-249 g and were housed >= 1 week before surgery)",
    weight_range   = "200-249 g body weight on arrival; weight on the experimental day was tested as a covariate in the screen and not retained in the final model",
    sex_female_pct = 0,
    race_ethnicity = NA,
    disease_state  = paste(
      "Two arms: naive saline-control adult Sprague-Dawley rats with",
      "a hippocampal microdialysis probe and femoral artery / vein",
      "cannulae; and kainate-induced status epilepticus rats studied",
      "7 days post-SE (a temporal-lobe-epilepsy paradigm in which",
      "P-gp expression at the BBB and on intra-parenchymal cells is",
      "hypothesised to be upregulated)."
    ),
    dose_range     = paste(
      "Three IV quinidine infusion arms in the control cohort:",
      "12.5 ug/min over 4 h (3000 ug cumulative), 25 ug/min over 4 h",
      "(6000 ug cumulative), and 100 ug/min over 30 min (3000 ug",
      "cumulative). Kainate cohort received the 100 ug/min over 30 min",
      "arm only. Half of each arm received tariquidar 15 mg/kg IP",
      "30 min before quinidine; the other half received 5% glucose /",
      "saline vehicle."
    ),
    regions        = "preclinical (in-vivo rat); Leiden University, The Netherlands",
    notes          = paste(
      "Of 74 rats initially enrolled, 10 were excluded for technical",
      "reasons (cannula clogging, microdialysis flow irregularities,",
      "death, or insufficient recovery after kainate) leaving 64",
      "evaluable. The PK-fit population is 58 of those 64 (the 6",
      "retrodialysis-only rats contributed in-vivo probe-recovery",
      "calibration data, not PK data). Table I (page 88) details the",
      "per-arm counts: control + vehicle 9 + 8 + 7 = 24 rats across",
      "the three quinidine arms; control + TQD 7 + 8 + 6 = 21 rats;",
      "kainate + vehicle 8 rats (100 ug/min over 30 min only);",
      "kainate + TQD 8 rats (100 ug/min over 30 min only). Total",
      "kainate-treated 16 of 64 (= 25%)."
    )
  )

  ini({
    # --------------------------------------------------------------
    # Plasma 2-compartment structural parameters. Values from
    # Syvanen 2012 Table II, "Estimation - 2898" column (final
    # model with brain compartments; the "Estimation plasma only - 695"
    # column is the intermediate plasma-only fit per the paper's
    # three-step Data Analysis procedure). Reference (typical)
    # values correspond to DIS_POSTSE_KAINATE = 0 and
    # CONMED_TARIQUIDAR = 0 (control rat + vehicle co-administration).
    lcl <- log(32.1);  label("Systemic plasma clearance CL (mL/min) -- reference: control rat + vehicle")   # Table II: CL = 32.1 mL/min (RSE 5.1%)
    lvc <- log(146);   label("Central plasma volume V1 (mL)")                                                # Table II: V1 = 146 mL (RSE 42%)
    lq  <- log(66.6);  label("Inter-compartmental plasma clearance Q (mL/min)")                              # Table II: Q2 = 66.6 mL/min (RSE 11%)
    lvp <- log(1840);  label("Peripheral plasma volume V2 (mL) -- reference: control rat")                  # Table II: V2 = 1840 mL (RSE 6.5%)

    # --------------------------------------------------------------
    # Brain transport structural parameters (paper-specific names
    # following the Xie_2000_m3g_rat.R brain-microdialysis precedent).
    # Q_out is FIXED per paper Results 'Development of the
    # Pharmacokinetic Model' (page 90): "Q_out was fixed to a value
    # of 10.8 mL/min" because in the bootstrap analysis V_Br1 was
    # occasionally estimated to an unrealistically high value while
    # the ratio V_Br1 / Q_out remained stable; the 10.8 mL/min anchor
    # is the value obtained when Q_out was left free in the final
    # model fit.
    lqout <- fixed(log(10.8));  label("BBB efflux clearance Q_out (mL/min), FIXED -- reference: vehicle co-administration")  # Table II: Q_out = 10.8 mL/min FIXED
    lf1   <- log(0.0619);       label("Ratio Q_in / Q_out (unitless)")                                                       # Table II: f1 = 0.0619 (RSE 8.2%)
    lvbr  <- log(269);          label("Brain ECF volume of distribution V_Br (mL) -- reference: control rat")                # Table II: V_Br = 269 mL (RSE 14%)
    lf2   <- log(8.67);         label("Ratio Cbrain_deep / Cbrain_csf (unitless) -- reference: vehicle co-administration")   # Table II: f2 = 8.67 (RSE 2.7%)

    # --------------------------------------------------------------
    # Covariate effects. The paper's NONMEM parameterisation is
    # multiplicative on the linear-scale typical-value parameters
    # (Table IV equations):
    #     theta_individual = theta_pop * theta_covar^COVARIATE * exp(eta)
    # In nlmixr2's additive log-scale form this becomes:
    #     param = exp(l<param> + e_<cov>_<param> * COVARIATE)
    # with e_<cov>_<param> = log(theta_covar). RSE quoted in the
    # source-trace comments is the paper-reported relative SE of
    # theta_covar (not of e_<cov>_<param>).
    e_dis_postse_kainate_cl      <- log(0.537);  label("Effect of DIS_POSTSE_KAINATE on log-CL (kainate decreases CL 1.86-fold)")       # Table II: theta_kainate(CL)        = 0.537 (RSE 8.5%)
    e_dis_postse_kainate_vp      <- log(0.678);  label("Effect of DIS_POSTSE_KAINATE on log-V2 (kainate decreases V2 1.47-fold)")       # Table II: theta_kainate(V2)        = 0.678 (RSE 7.0%)
    e_dis_postse_kainate_vbr     <- log(0.592);  label("Effect of DIS_POSTSE_KAINATE on log-V_Br (kainate decreases V_Br 1.69-fold)")   # Table II: theta_kainate(V_Br)      = 0.592 (RSE 13%)
    e_conmed_tariquidar_cl   <- log(0.835);  label("Effect of CONMED_TARIQUIDAR on log-CL (TQD decreases CL 1.20-fold)")        # Table II: theta_tariquidar(CL)      = 0.835 (RSE 6.8%)
    e_conmed_tariquidar_qout <- log(0.366);  label("Effect of CONMED_TARIQUIDAR on log-Q_out (TQD decreases Q_out 2.73-fold)")  # Table II: theta_tariquidar(Q_out)   = 0.366 (RSE 14%)
    e_conmed_tariquidar_qin  <- log(2.65);   label("Effect of CONMED_TARIQUIDAR on log-Q_in (TQD increases Q_in 2.65-fold)")    # Table II: theta_tariquidar(Q_in)    = 2.65  (RSE 17%)
    e_conmed_tariquidar_f2   <- log(5.66);   label("Effect of CONMED_TARIQUIDAR on log-f2 (TQD increases f2 5.66-fold)")        # Table II: theta_tariquidar(f2)      = 5.66  (RSE 33%)

    # --------------------------------------------------------------
    # Inter-animal variability. Syvanen 2012 Methods 'Data Analysis'
    # Eq. 3 specifies the exponential variance model:
    #     theta_i = theta_pop * exp(eta_i),  eta_i ~ N(0, omega^2)
    # Table II reports the omega^2 values directly under the
    # eta(parameter) labels (NONMEM $OMEGA convention); they are
    # entered here as variances on the log scale. The very large
    # eta(V1) = 1.92 variance (Table II shrinkage 18%, bootstrap 2.03)
    # reflects high inter-animal uncertainty on central plasma volume
    # in this design (rich brain-ECF sampling but only sparse early-
    # phase arterial blood sampling that constrains V1); the variance
    # is taken at face value from the paper.
    etalcl  ~ 0.0682   # Table II: omega^2(CL)   = 0.0682 (RSE 20%, shrinkage  6%); approximate CV 27%
    etalvc  ~ 1.92     # Table II: omega^2(V1)   = 1.92   (RSE 30%, shrinkage 18%); approximate CV 240% -- high IIV on central plasma volume per paper Table II
    etalvp  ~ 0.0686   # Table II: omega^2(V2)   = 0.0686 (RSE 34%, shrinkage 15%); approximate CV 27%
    etalqin ~ 0.101    # Table II: omega^2(Q_in) = 0.101  (RSE 19%, shrinkage 13%); approximate CV 33%
    etalvbr ~ 0.196    # Table II: omega^2(V_Br) = 0.196  (RSE 31%, shrinkage 16%); approximate CV 47%

    # --------------------------------------------------------------
    # Residual error: independent proportional models for each of
    # the three observation streams. Table II reports sigma^2 (the
    # variance on the multiplicative noise term in the model
    # y_obs = y_pred * (1 + eps), eps ~ N(0, sigma^2)). propSd is
    # entered as sqrt(sigma^2) per nlmixr2's convention that
    # ~prop(propSd) takes the SD.
    propSd             <- sqrt(0.270);  label("Proportional residual SD on plasma Cc (fraction)")                       # Table II: sigma^2(plasma)      = 0.270 (RSE  7.0%) -> SD 0.520
    propSd_Cbrain_csf  <- sqrt(0.185);  label("Proportional residual SD on brain ECF Cbrain_csf (fraction)")            # Table II: sigma^2(brain ECF)   = 0.185 (RSE  9.4%) -> SD 0.430
    propSd_Cbrain_deep <- sqrt(0.694);  label("Proportional residual SD on total-brain Cbrain_deep (fraction)")         # Table II: sigma^2(total brain) = 0.694 (RSE 21%)   -> SD 0.833
  })

  model({
    # --------------------------------------------------------------
    # Individual structural parameters with covariate effects. The
    # exponential transform expands the paper's multiplicative form
    # (Table IV equations) into the additive log-scale form used
    # by nlmixr2.
    cl   <- exp(lcl   + e_dis_postse_kainate_cl      * DIS_POSTSE_KAINATE   + e_conmed_tariquidar_cl   * CONMED_TARIQUIDAR + etalcl)
    vc   <- exp(lvc   + etalvc)
    q    <- exp(lq)
    vp   <- exp(lvp   + e_dis_postse_kainate_vp      * DIS_POSTSE_KAINATE   + etalvp)
    vbr  <- exp(lvbr  + e_dis_postse_kainate_vbr     * DIS_POSTSE_KAINATE   + etalvbr)
    qout <- exp(lqout + e_conmed_tariquidar_qout * CONMED_TARIQUIDAR)
    f1   <- exp(lf1)
    # Table IV: Q_in = f1 * THETA8 * THETA13^TARIQUIDAR * exp(eta5)
    # where THETA8 is the BASE Q_out value (10.8 mL/min, FIXED)
    # BEFORE its own tariquidar modification THETA9. The equivalent
    # additive log-scale construction passes the unmodified lqout
    # into qin and adds the qin-specific tariquidar effect.
    qin  <- f1 * exp(lqout + e_conmed_tariquidar_qin * CONMED_TARIQUIDAR + etalqin)
    f2   <- exp(lf2 + e_conmed_tariquidar_f2 * CONMED_TARIQUIDAR)

    # --------------------------------------------------------------
    # Concentrations driving the ODE and observation equations.
    # Compartments hold amounts in ng; volumes in mL; concentrations
    # therefore arrive in ng/mL natively (matching the paper's
    # reported units without an explicit conversion factor).
    Cc          <- central     / vc
    Cp          <- peripheral1 / vp
    Cbrain_csf  <- brain_csf   / vbr
    Cbrain_deep <- f2 * Cbrain_csf

    # --------------------------------------------------------------
    # ODE system (Syvanen 2012 Fig. 2 schematic). Quinidine is dosed
    # by IV infusion into the central compartment (the events table
    # carries cmt = central, evid = 1, rate / amt records). The brain
    # ECF compartment exchanges bidirectionally with central via
    # Q_in / Q_out; the deep brain compartment is NOT an ODE state
    # because the paper could not estimate its rate constants and
    # modelled it as the algebraic equilibrium ratio
    # Cbrain_deep = f2 * Cbrain_csf (paper Results 'Development of
    # the Pharmacokinetic Model' and Discussion p93: "even if brain
    # samples had been obtained at different time points it would
    # still have been possible to use the present model structure").
    d/dt(central)     <- -cl * Cc - q * Cc + q * Cp - qin * Cc + qout * Cbrain_csf
    d/dt(peripheral1) <-  q * Cc - q * Cp
    d/dt(brain_csf)   <-  qin * Cc - qout * Cbrain_csf

    # --------------------------------------------------------------
    # Observation equations: three independent proportional-error
    # streams matching Table II's three sigma^2 values.
    Cc          ~ prop(propSd)
    Cbrain_csf  ~ prop(propSd_Cbrain_csf)
    Cbrain_deep ~ prop(propSd_Cbrain_deep)
  })
}
