Kervezee_2014_quinidine_rat <- function() {
  description <- paste(
    "PBPK (semi-mechanistic brain-distribution, unbound-only). Preclinical",
    "(rat). Nine-compartment brain-distribution model for i.v. quinidine in",
    "male Wistar rats (Kervezee 2014), fitted to unbound plasma, brain",
    "extracellular fluid (ECF), cerebrospinal fluid at the lateral ventricle",
    "(CSF-LV) and cisterna magna (CSF-CM), and total deep-brain tissue.",
    "Topology from Westerhout 2013: plasma (central), two peripheral tissue",
    "compartments (peripheral1, peripheral2), one deep-brain intracellular",
    "compartment (brain_deep), one brain ECF compartment (brain_ecf), and",
    "four sequential CSF sub-compartments (csf_lv, csf_tfv, csf_cm, csf_sas)",
    "draining back into plasma via CSF bulk flow. P-glycoprotein-mediated",
    "transport enters as an additional clearance component at the BBB",
    "(subtracted from PL-to-brain influx, added to brain-to-PL efflux) and",
    "on plasma elimination. Kervezee 2014 introduces the diurnal-period",
    "covariate PERIOD_ACTIVE (0 = resting / lights-on, 1 = active / lights-",
    "off) that acts on five parameters: the P-gp components of CL_DBR-PL,",
    "CL_PL-ECF, CL_ECF-PL, and CL_PL-LV, plus CSF bulk flow Q_CSF. Brain",
    "compartment volumes and Q_ECF are fixed to physiological values from",
    "Westerhout 2013 refs 38-47; plasma volume V_PL is fixed to the rat",
    "plasma volume; peripheral volumes are estimated.")
  reference <- paste(
    "Kervezee L, Hartman R, van den Berg DJ, Shimizu S, Emoto-Yamamoto Y,",
    "Meijer JH, de Lange ECM. Diurnal variation in P-glycoprotein-mediated",
    "transport and cerebrospinal fluid turnover in the brain. AAPS J.",
    "2014;16(5):1029-1037. doi:10.1208/s12248-014-9625-4.",
    "PBPK topology and rate-constant conventions inherited from Westerhout",
    "J, Smeets J, Danhof M, de Lange ECM. The impact of P-gp functionality",
    "on non-steady state relationships between CSF and brain extracellular",
    "fluid. J Pharmacokinet Pharmacodyn. 2013;40:327-342.",
    "doi:10.1007/s10928-013-9314-4.")
  vignette <- "Kervezee_2014_quinidine_rat"
  units <- list(
    time          = "min",
    dosing        = "mg",
    concentration = "ng/mL"
  )

  paper_specific_compartments <- c(
    "brain_deep", "brain_ecf", "csf_lv", "csf_tfv", "csf_cm", "csf_sas"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. analyte/specimen proposed by a local model from the
  # model description; units derived from the units block. verified = FALSE
  # means NOT checked against the source paper.
  compartmentData <- list(
    central     = list(analyte = "quinidine unbound", units = "mg", specimen = "plasma", verified = FALSE),
    peripheral1 = list(analyte = "quinidine unbound", units = "mg", specimen = "plasma", verified = FALSE),
    peripheral2 = list(analyte = "quinidine unbound", units = "mg", specimen = "plasma", verified = FALSE),
    brain_deep  = list(analyte = "quinidine unbound", units = "mg", specimen = "tissue", verified = FALSE),
    brain_ecf   = list(analyte = "quinidine unbound", units = "mg", specimen = "tissue", verified = FALSE),
    csf_lv      = list(analyte = "quinidine unbound", units = "mg", specimen = "CSF", verified = FALSE),
    csf_tfv     = list(analyte = "quinidine unbound", units = "mg", specimen = "CSF", verified = FALSE),
    csf_cm      = list(analyte = "quinidine unbound", units = "mg", specimen = "CSF", verified = FALSE),
    csf_sas     = list(analyte = "quinidine unbound", units = "mg", specimen = "CSF", verified = FALSE)
  )

  covariateData <- list(
    PERIOD_ACTIVE = list(
      description        = paste(
        "Diurnal-period indicator: 1 = active period (lights-off / dark",
        "phase, ZT12-ZT24) for the nocturnally-active Wistar rat; 0 =",
        "resting period (lights-on / light phase, ZT0-ZT12). Selects the",
        "active-period vs resting-period estimates of the P-gp component",
        "of four BBB clearances (CL_DBR-PL, CL_PL-ECF, CL_ECF-PL, CL_PL-",
        "LV) and of CSF bulk flow Q_CSF."),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (resting period; lights-on)",
      notes              = paste(
        "Time-fixed per animal in the source study (each rat is dosed",
        "at a single ZT). Kervezee 2014 aggregated six experimental ZT",
        "levels (ZT0, 4, 8 as resting; ZT12, 16, 20 as active) into the",
        "two-level PERIOD_ACTIVE variable for the final covariate model",
        "(Table I). For nocturnally-active species (rats, mice)",
        "PERIOD_ACTIVE = 1 corresponds to the dark / behaviourally-",
        "active phase; for diurnally-active species (humans) the",
        "mapping is inverted."),
      source_name        = "time of administration"
    )
  )

  population <- list(
    species        = "rat (Wistar, male)",
    n_subjects     = 55L,
    n_studies      = 1L,
    age_range      = "Adult; specific age not reported",
    weight_range   = "Not reported (source paper Methods: 'Male Wistar rats (Charles River, The Netherlands)')",
    sex_female_pct = 0,
    race_ethnicity = NA,
    disease_state  = paste(
      "Healthy male Wistar rats (Charles River, The Netherlands) housed",
      "under standard 12:12 light-dark cycle with free access to food",
      "and water. Two sub-experiments contribute to the fitted dataset:",
      "(1) brain distribution (N = 40; n=6-8 per ZT level; blood + total-",
      "brain sampling at t = -10, 10, 30, 60, 90, 120, 150, 180, 210,",
      "240 min); (2) intracerebral microdialysis (N = 15; ZT8 and ZT20;",
      "microdialysis probes in caudate putamen for ECF and cisterna",
      "magna for CSF-CM, plus lateral-ventricle CSF sampling from",
      "Westerhout 2013 pooled data). Microdialysis probe recovery from",
      "the retrodialysis calibration was 13 +/- 1.4% (4 mm CP probe)",
      "and 8.4 +/- 2.6% (1 mm CM probe)."
    ),
    dose_range     = paste(
      "10 mg/kg i.v. infusion over 10 min at t = 0 (rate 250 uL/min/kg",
      "in 5% glucose vehicle). Pre-treatment either vehicle (5% glucose",
      "in saline, N ~= 40) or 15 mg/kg i.v. tariquidar (N ~= 15) at t =",
      "-25 min over 10 min at 500 uL/min/kg (Methods, Drug",
      "Administration section). Tariquidar arms abolish the P-gp",
      "component of brain clearance and are not used to identify the",
      "diurnal-period effect; they contribute to the pooled likelihood",
      "for identifiability of the passive-clearance parameters."
    ),
    regions        = "Leiden, the Netherlands (single laboratory; DEC12088)",
    notes          = paste(
      "Demographics from Kervezee 2014 Materials and Methods (Animals;",
      "Drug Administration, Serial Blood Sampling, and Collection of",
      "Brain Tissue; Intracerebral Microdialysis; Supplemental Table",
      "I). N=55 is the pooled Kervezee 2014 vehicle + tariquidar cohort",
      "reported in Supplemental Table I; the PBPK fit additionally",
      "pooled equivalent-method microdialysis data from Westerhout et",
      "al. 2013 (a preceding non-diurnal study using the same probes",
      "and dosing regimen in the resting period) for identifiability",
      "of the passive-clearance parameters. Total-brain concentrations",
      "are measured in brain homogenate at t = 240 min after",
      "transcardial PBS perfusion; unbound plasma concentrations are",
      "derived from total plasma corrected by the unbound fraction",
      "(fu = 0.286 +/- 0.006, not model-parameterised); microdialysate",
      "concentrations are corrected to free-drug ECF / CSF",
      "concentrations by probe-specific in vivo recovery."
    )
  )

  ini({
    # ----------------------------------------------------------------------
    # Structural parameters -- Kervezee 2014 Table II (page 1035).
    # ODE topology, rate-constant conventions, and the P-gp sign-of-effect
    # for BBB compartments are inherited from Westerhout 2013 Appendix
    # (JPKPD 40:327-342, differential equations page 339-340):
    #   k_PL->X = (CL_PL-X,p  -  CL_PL-X,Pgp) / V_PL      [P-gp reduces influx]
    #   k_X->PL = (CL_X-PL,p  +  CL_X-PL,Pgp) / V_X       [P-gp increases efflux]
    #   k_E     = (CL_E,p     +  CL_E,Pgp)    / V_PL      [elimination sum]
    # Applied to X in {DBR (deep brain), ECF, LV (lateral ventricle),
    # TFV (third+fourth ventricle), CM (cisterna magna), SAS (subarachnoid
    # space)}. CSF bulk flow follows the sequential chain
    # ECF -> LV -> TFV -> CM -> SAS -> PL at rates (Q_ECF / V_ECF) and
    # (Q_CSF / V_upstream).
    #
    # Units: all clearances converted to mL/min, all volumes to mL, all rate
    # constants therefore 1/min. Paper's mixed uL/min / mL/min / L / uL
    # values are preserved in the source-trace comments for auditability.
    # Non-canonical parameter names (lclPLDBR_p, lclDBRPL_gr, ...) reflect
    # the paper's per-transition passive-vs-Pgp / resting-vs-active split;
    # justified in the vignette Assumptions and deviations.
    # ----------------------------------------------------------------------

    # --- Plasma elimination -------------------------------------------------
    lclE_p <- log(83.3);        label("Passive plasma elimination CL_E,p (mL/min)")      # Kervezee 2014 Table II: CL_E Passive 83.3 mL/min, RSE 10.3%
    lclE_g <- log(52.7);        label("P-gp plasma elimination CL_E,Pgp (mL/min)")       # Kervezee 2014 Table II: CL_E +P-gp 52.7 mL/min, RSE 6.5%

    # --- Peripheral distribution -------------------------------------------
    lqPer1 <- log(831);         label("Plasma-peripheral1 flow Q_PL-PER1 (mL/min)")      # Kervezee 2014 Table II: Q_PL-PER1 Passive 831 mL/min, RSE 9.6%
    lqPer2 <- log(93.5);        label("Plasma-peripheral2 flow Q_PL-PER2 (mL/min)")      # Kervezee 2014 Table II: Q_PL-PER2 Passive 93.5 mL/min, RSE 18.3%

    # --- Deep brain (DBR) -- BBB transport ---------------------------------
    lclPLDBR_p  <- log(0.982);  label("Passive PL->DBR clearance (mL/min)")              # Kervezee 2014 Table II: CL_PL-Deep Brain Passive 982 uL/min = 0.982 mL/min, RSE 14.7%; +P-gp NA (fixed to 0 below)
    lclDBRPL_p  <- log(0.0123); label("Passive DBR->PL clearance (mL/min)")              # Kervezee 2014 Table II: CL_Deep Brain-PL Passive 12.3 uL/min = 0.0123 mL/min, RSE 19.8%
    lclDBRPL_gr <- log(0.228);  label("P-gp DBR->PL clearance, resting period (mL/min)") # Kervezee 2014 Table II: CL_Deep Brain-PL +P-gp Resting 228 uL/min = 0.228 mL/min, RSE 17.7% (PERIOD_ACTIVE = 0)
    lclDBRPL_ga <- log(0.659);  label("P-gp DBR->PL clearance, active period (mL/min)")  # Kervezee 2014 Table II: CL_Deep Brain-PL +P-gp Active 659 uL/min = 0.659 mL/min, RSE 15.6% (PERIOD_ACTIVE = 1)

    # --- Brain ECF -- BBB transport ----------------------------------------
    lclPLECF_p  <- log(0.0257); label("Passive PL->ECF clearance (mL/min)")              # Kervezee 2014 Table II: CL_PL-ECF Passive 25.7 uL/min = 0.0257 mL/min, RSE 12.5%
    lclPLECF_gr <- log(0.017);  label("P-gp PL->ECF clearance, resting period (mL/min)") # Kervezee 2014 Table II: CL_PL-ECF +P-gp Resting 17 uL/min = 0.017 mL/min, RSE 18.5%
    lclPLECF_ga <- log(0.0183); label("P-gp PL->ECF clearance, active period (mL/min)")  # Kervezee 2014 Table II: CL_PL-ECF +P-gp Active 18.3 uL/min = 0.0183 mL/min, RSE 18.5%
    lclECFPL_p  <- log(0.00463); label("Passive ECF->PL clearance (mL/min)")             # Kervezee 2014 Table II: CL_ECF-PL Passive 4.63 uL/min = 0.00463 mL/min, RSE 15.2%
    lclECFPL_gr <- log(0.00314); label("P-gp ECF->PL clearance, resting period (mL/min)")# Kervezee 2014 Table II: CL_ECF-PL +P-gp Resting 3.14 uL/min = 0.00314 mL/min, RSE 36.3%
    lclECFPL_ga <- log(0.00398); label("P-gp ECF->PL clearance, active period (mL/min)") # Kervezee 2014 Table II: CL_ECF-PL +P-gp Active 3.98 uL/min = 0.00398 mL/min, RSE 48.0%

    # --- CSF lateral ventricle (LV) -- BCSFB transport ---------------------
    lclPLLV_p  <- log(0.00323); label("Passive PL->CSF-LV clearance (mL/min)")           # Kervezee 2014 Table II: CL_PL-LV Passive 3.23 uL/min = 0.00323 mL/min, RSE 19.9%
    lclPLLV_gr <- log(0.00155); label("P-gp PL->CSF-LV clearance, resting (mL/min)")     # Kervezee 2014 Table II: CL_PL-LV +P-gp Resting 1.55 uL/min = 0.00155 mL/min, RSE 30.1%
    lclPLLV_ga <- log(0.00244); label("P-gp PL->CSF-LV clearance, active (mL/min)")      # Kervezee 2014 Table II: CL_PL-LV +P-gp Active 2.44 uL/min = 0.00244 mL/min, RSE 17.5%
    lclLVPL_p  <- log(0.000513); label("Passive CSF-LV->PL clearance (mL/min)")          # Kervezee 2014 Table II: CL_LV-PL Passive 0.513 uL/min = 0.000513 mL/min, RSE 24.0%; +P-gp NA (fixed to 0 below)

    # --- CSF cisterna magna (CM) -- BCSFB transport (no P-gp) --------------
    lclPLCM_p  <- log(0.000753); label("Passive PL->CSF-CM clearance (mL/min)")          # Kervezee 2014 Table II: CL_PL-CM Passive 0.753 uL/min = 0.000753 mL/min, RSE 23.5%; +P-gp NA (fixed to 0 below)
    lclCMPL_p  <- log(0.00102);  label("Passive CSF-CM->PL clearance (mL/min)")          # Kervezee 2014 Table II: CL_CM-PL Passive 1.02 uL/min = 0.00102 mL/min, RSE 33.7%; +P-gp NA (fixed to 0 below)

    # --- ECF / CSF bulk flows ----------------------------------------------
    lqEcf   <- fixed(log(0.0002));  label("Brain ECF bulk flow Q_ECF (mL/min)")          # Kervezee 2014 Table II: Q_ECF 0.2 uL/min = 0.0002 mL/min, physiological (Westerhout 2013 refs 39, 46); fixed
    lqCsf_r <- log(0.000522);       label("CSF bulk flow Q_CSF, resting (mL/min)")       # Kervezee 2014 Table II: Q_CSF Resting 0.522 uL/min = 0.000522 mL/min, RSE 28.5%
    lqCsf_a <- log(0.000227);       label("CSF bulk flow Q_CSF, active (mL/min)")        # Kervezee 2014 Table II: Q_CSF Active 0.227 uL/min = 0.000227 mL/min, RSE 36.0%

    # --- Compartment volumes -----------------------------------------------
    # Physiological volumes fixed to rat anatomy per Westerhout 2013 Methods
    # (page 331) and Kervezee 2014 Table II 'b' superscript (physiological
    # values from literature).
    lvPl   <- fixed(log(10.6));    label("Plasma volume V_PL (mL)")                       # Kervezee 2014 Table II: V_PL 10.6 mL (physiological, fixed)
    lvPer1 <- log(7420);            label("Peripheral1 volume V_PER1 (mL)")               # Kervezee 2014 Table II: V_PER1 7.42 L = 7420 mL, RSE 5.7%
    lvPer2 <- log(7090);            label("Peripheral2 volume V_PER2 (mL)")               # Kervezee 2014 Table II: V_PER2 7.09 L = 7090 mL, RSE 17.3%
    lvDbr  <- fixed(log(1.44));    label("Deep brain volume V_DBR (mL)")                  # Kervezee 2014 Table II: V_Deep Brain 1440 uL = 1.44 mL (physiological, fixed; Westerhout 2013 ref 38)
    lvEcf  <- fixed(log(0.29));    label("Brain ECF volume V_ECF (mL)")                   # Kervezee 2014 Table II: V_ECF 290 uL = 0.29 mL (physiological, fixed; Westerhout 2013 ref 39)
    lvLv   <- fixed(log(0.05));    label("CSF lateral-ventricle volume V_LV (mL)")        # Kervezee 2014 Table II: V_LV 50 uL = 0.05 mL (physiological, fixed; Westerhout 2013 refs 41, 42)
    lvTfv  <- fixed(log(0.05));    label("CSF third+fourth-ventricle volume V_TFV (mL)")  # Kervezee 2014 Table II: V_TFV 50 uL = 0.05 mL (physiological, fixed; Westerhout 2013 ref 43)
    lvCm   <- fixed(log(0.017));   label("CSF cisterna-magna volume V_CM (mL)")           # Kervezee 2014 Table II: V_CM 17 uL = 0.017 mL (physiological, fixed; Westerhout 2013 refs 44, 45)
    lvSas  <- fixed(log(0.18));    label("CSF subarachnoid-space volume V_SAS (mL)")      # Kervezee 2014 Table II: V_SAS 180 uL = 0.18 mL (physiological, fixed; Westerhout 2013 refs 40, 43)

    # ----------------------------------------------------------------------
    # Inter-individual variability. Kervezee 2014 Table II reports IIV
    # only on the plasma elimination CL_E; all other parameters are
    # estimated as typical values (no IIV) or fixed to physiology. The
    # reported CV% is interpreted as the log-normal coefficient of
    # variation of the individual-value distribution:
    #   sigma^2 = log(1 + CV^2) = log(1 + 0.332^2) = log(1.1102) = 0.1046
    # The eta multiplies the total (passive + P-gp) elimination clearance
    # via cl_e_typ = exp(lclE_p) + exp(lclE_g); cl_e = cl_e_typ *
    # exp(etalclE), preserving the Table II split for source-tracing.
    # ----------------------------------------------------------------------
    etalclE ~ 0.1046                                                                       # Kervezee 2014 Table II: IIV CL_E 33.2% CV, RSE 17.2% -> sigma^2 = log(1 + 0.332^2) = 0.1046

    # ----------------------------------------------------------------------
    # Residual error. Kervezee 2014 Table II reports proportional
    # residual CV% for the five observation types (PL, ECF, LV, CM,
    # Deep brain). Encoded as bare proportional-SD parameters (CV / 100
    # gives the propSd in nlmixr2's proportional form since concentration
    # ~ prop(sd) implements Y = Y_pred * (1 + eps * sd)).
    # ----------------------------------------------------------------------
    propSd            <- 0.428;   label("Proportional residual SD, unbound plasma (fraction)")  # Kervezee 2014 Table II: Residual error PL 42.8% CV, RSE 13.9%
    propSd_Cbrain_ecf  <- 0.330;   label("Proportional residual SD, brain ECF (fraction)")       # Kervezee 2014 Table II: Residual error ECF 33.0% CV, RSE 13.9%
    propSd_Ccsf_lv     <- 0.319;   label("Proportional residual SD, CSF LV (fraction)")          # Kervezee 2014 Table II: Residual error LV 31.9% CV, RSE 18.7%
    propSd_Ccsf_cm     <- 0.362;   label("Proportional residual SD, CSF CM (fraction)")          # Kervezee 2014 Table II: Residual error CM 36.2% CV, RSE 13.8%
    propSd_Cbrain_deep <- 0.356;   label("Proportional residual SD, deep brain (fraction)")      # Kervezee 2014 Table II: Residual error Deep brain 35.6% CV, RSE 13.4%
  })

  model({
    # ------------------------------------------------------------------
    # 1. Convert log-transformed parameters into linear rate/volume space.
    #    Individual plasma elimination applies the CL_E IIV to the
    #    passive + P-gp sum (Table II reports a single IIV entry for
    #    CL_E covering the combined elimination).
    # ------------------------------------------------------------------
    cl_e_typ <- exp(lclE_p) + exp(lclE_g)
    cl_e     <- cl_e_typ * exp(etalclE)

    q_per1   <- exp(lqPer1)
    q_per2   <- exp(lqPer2)

    # PL <-> DBR: passive both directions; +P-gp increases efflux only
    # (paper reports NA for the influx +P-gp entry -> fixed to zero).
    cl_pl_dbr_p  <- exp(lclPLDBR_p)
    cl_dbr_pl_p  <- exp(lclDBRPL_p)
    cl_dbr_pl_gr <- exp(lclDBRPL_gr)
    cl_dbr_pl_ga <- exp(lclDBRPL_ga)
    cl_dbr_pl_g  <- cl_dbr_pl_gr * (1 - PERIOD_ACTIVE) +
                    cl_dbr_pl_ga *      PERIOD_ACTIVE

    # PL <-> ECF: passive both directions; +P-gp in both directions
    # (subtracts on influx, adds on efflux per Westerhout 2013).
    cl_pl_ecf_p  <- exp(lclPLECF_p)
    cl_pl_ecf_gr <- exp(lclPLECF_gr)
    cl_pl_ecf_ga <- exp(lclPLECF_ga)
    cl_pl_ecf_g  <- cl_pl_ecf_gr * (1 - PERIOD_ACTIVE) +
                    cl_pl_ecf_ga *      PERIOD_ACTIVE
    cl_ecf_pl_p  <- exp(lclECFPL_p)
    cl_ecf_pl_gr <- exp(lclECFPL_gr)
    cl_ecf_pl_ga <- exp(lclECFPL_ga)
    cl_ecf_pl_g  <- cl_ecf_pl_gr * (1 - PERIOD_ACTIVE) +
                    cl_ecf_pl_ga *      PERIOD_ACTIVE

    # PL <-> LV: passive both directions; +P-gp on influx only
    # (paper reports NA for the LV->PL +P-gp entry -> fixed to zero).
    cl_pl_lv_p  <- exp(lclPLLV_p)
    cl_pl_lv_gr <- exp(lclPLLV_gr)
    cl_pl_lv_ga <- exp(lclPLLV_ga)
    cl_pl_lv_g  <- cl_pl_lv_gr * (1 - PERIOD_ACTIVE) +
                   cl_pl_lv_ga *      PERIOD_ACTIVE
    cl_lv_pl_p  <- exp(lclLVPL_p)

    # PL <-> CM: passive-only both directions (paper reports NA for
    # both +P-gp entries -> fixed to zero, no diurnal effect).
    cl_pl_cm_p <- exp(lclPLCM_p)
    cl_cm_pl_p <- exp(lclCMPL_p)

    # ECF and CSF bulk flows.
    q_ecf   <- exp(lqEcf)
    q_csf_r <- exp(lqCsf_r)
    q_csf_a <- exp(lqCsf_a)
    q_csf   <- q_csf_r * (1 - PERIOD_ACTIVE) + q_csf_a * PERIOD_ACTIVE

    # Volumes.
    vPl   <- exp(lvPl)
    vPer1 <- exp(lvPer1)
    vPer2 <- exp(lvPer2)
    vDbr  <- exp(lvDbr)
    vEcf  <- exp(lvEcf)
    vLv   <- exp(lvLv)
    vTfv  <- exp(lvTfv)
    vCm   <- exp(lvCm)
    vSas  <- exp(lvSas)

    # ------------------------------------------------------------------
    # 2. Effective rate constants -- Westerhout 2013 Appendix.
    #    P-gp subtracts on influx and adds on efflux.
    # ------------------------------------------------------------------
    kE      <- cl_e / vPl

    kPlPer1 <- q_per1 / vPl
    kPer1Pl <- q_per1 / vPer1
    kPlPer2 <- q_per2 / vPl
    kPer2Pl <- q_per2 / vPer2

    kPlDbr  <-  cl_pl_dbr_p                 / vPl
    kDbrPl  <- (cl_dbr_pl_p + cl_dbr_pl_g)  / vDbr

    kPlEcf  <- (cl_pl_ecf_p - cl_pl_ecf_g)  / vPl
    kEcfPl  <- (cl_ecf_pl_p + cl_ecf_pl_g)  / vEcf

    kPlLv   <- (cl_pl_lv_p  - cl_pl_lv_g)   / vPl
    kLvPl   <-  cl_lv_pl_p                  / vLv

    kPlCm   <-  cl_pl_cm_p                  / vPl
    kCmPl   <-  cl_cm_pl_p                  / vCm

    kEcfLv  <- q_ecf / vEcf
    kLvTfv  <- q_csf / vLv
    kTfvCm  <- q_csf / vTfv
    kCmSas  <- q_csf / vCm
    kSasPl  <- q_csf / vSas

    # ------------------------------------------------------------------
    # 3. ODE system -- Westerhout 2013 Appendix (mass-balance form):
    #      Central = plasma unbound pool (i.v. dose lands here).
    #      Peripheral1, Peripheral2 = passive tissue distribution.
    #      brain_deep = deep brain intracellular + tissue pool.
    #      brain_ecf  = brain extracellular fluid.
    #      csf_lv -> csf_tfv -> csf_cm -> csf_sas -> central via Q_CSF.
    #    Note: brain_ecf -> csf_lv coupling via Q_ECF; brain_deep
    #    exchanges only with plasma (no direct brain-brain paths).
    # ------------------------------------------------------------------
    d/dt(central) <-
      (kPer1Pl * peripheral1) + (kPer2Pl * peripheral2) +
      (kDbrPl * brain_deep) + (kEcfPl * brain_ecf) +
      (kLvPl * csf_lv) + (kCmPl * csf_cm) + (kSasPl * csf_sas) -
      (kE + kPlPer1 + kPlPer2 + kPlDbr + kPlEcf + kPlLv + kPlCm) * central

    d/dt(peripheral1) <- kPlPer1 * central - kPer1Pl * peripheral1
    d/dt(peripheral2) <- kPlPer2 * central - kPer2Pl * peripheral2

    d/dt(brain_deep)  <- kPlDbr * central - kDbrPl * brain_deep

    d/dt(brain_ecf)   <- kPlEcf * central - (kEcfPl + kEcfLv) * brain_ecf

    d/dt(csf_lv)  <- kPlLv * central + kEcfLv * brain_ecf - (kLvPl + kLvTfv) * csf_lv
    d/dt(csf_tfv) <- kLvTfv * csf_lv - kTfvCm * csf_tfv
    d/dt(csf_cm)  <- kPlCm * central + kTfvCm * csf_tfv - (kCmPl + kCmSas) * csf_cm
    d/dt(csf_sas) <- kCmSas * csf_cm - kSasPl * csf_sas

    # ------------------------------------------------------------------
    # 4. Observations -- unbound concentrations in each compartment.
    #    Volumes are in mL, amounts in mg (i.v. dose in mg), so
    #    amount/volume = mg/mL = ug/uL = 1e6 ng/mL. Multiply by 1e6 to
    #    express concentrations in ng/mL (the paper's assay unit).
    # ------------------------------------------------------------------
    Cc          <- central     / vPl  * 1e6
    Cbrain_ecf  <- brain_ecf   / vEcf * 1e6
    Ccsf_lv     <- csf_lv      / vLv  * 1e6
    Ccsf_cm     <- csf_cm      / vCm  * 1e6
    Cbrain_deep <- brain_deep  / vDbr * 1e6

    Cc          ~ prop(propSd)
    Cbrain_ecf  ~ prop(propSd_Cbrain_ecf)
    Ccsf_lv     ~ prop(propSd_Ccsf_lv)
    Ccsf_cm     ~ prop(propSd_Ccsf_cm)
    Cbrain_deep ~ prop(propSd_Cbrain_deep)
  })
}
