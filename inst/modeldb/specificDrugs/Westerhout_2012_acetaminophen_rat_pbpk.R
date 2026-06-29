Westerhout_2012_acetaminophen_rat_pbpk <- function() {
  description <- paste(
    "PBPK (semi-mechanistic, regional brain) population PK model for",
    "acetaminophen (paracetamol) in 24 male Wistar WU rats (225-275 g),",
    "developed to investigate regional brain distribution kinetics with",
    "simultaneous microdialysis sampling in striatum (brain extracellular",
    "fluid), lateral ventricle (CSF_LV), and cisterna magna (CSF_CM) after",
    "a 10-min intravenous infusion of 15 mg/kg acetaminophen (Westerhout",
    "et al. 2012, AAPS J). Seven physiological compartments: plasma plus",
    "peripheral tissue plus brain extracellular fluid (with the brain",
    "intracellular space volume added per paper text page 5) plus four",
    "anatomically distinct CSF subcompartments (lateral ventricle, third",
    "and fourth ventricle combined, cisterna magna, subarachnoid space).",
    "Brain compartment volumes (V_pl, V_ICS, V_ECF, V_LV, V_TFV, V_CM,",
    "V_SAS) and bulk fluid flows (Q_ECF, Q_CSF) are fixed to literature",
    "physiological values for a 250-g rat; plasma-to-region and region-",
    "to-plasma clearances are estimated. The plasma-to-third-fourth-",
    "ventricle clearances CL15 / CL51 are structurally assumed equal to",
    "the plasma-to-lateral-ventricle clearances CL14 / CL41 (paper",
    "Results). An enterohepatic-recirculation continuous input",
    "F_abs * DOSE adds drug back to plasma to capture the apparent plateau",
    "after t = 120 min. The model is fitted to unbound plasma concentrations",
    "(plasma concentrations corrected to free fraction fu_p = 0.805); fu_p",
    "is documented in population metadata but does not enter the structural",
    "ODEs (see vignette Assumptions and deviations)."
  )
  reference <- paste(
    "Westerhout J, Ploeger B, Smeets J, Danhof M, de Lange ECM.",
    "Physiologically Based Pharmacokinetic Modeling to Investigate",
    "Regional Brain Distribution Kinetics in Rats. AAPS J.",
    "2012;14(3):543-553. doi:10.1208/s12248-012-9366-1.",
    sep = " "
  )
  vignette <- "Westerhout_2012_acetaminophen_rat_pbpk"

  # Five paper-mechanistic regional brain compartments: brain ECF (named
  # to match the paper notation; distinct from the canonical brain_csf
  # which is the deprecated alias of a generic brain CSF compartment) and
  # four anatomically distinct CSF subcompartments. None of these exist
  # as canonical brain_<region> entries in inst/references/compartment-
  # names.md, so they are declared here as paper-mechanistic.
  paper_specific_compartments <- c(
    "brain_ecf",
    "csf_lv", "csf_tfv", "csf_cm", "csf_sas"
  )

  units <- list(
    time          = "min",
    dosing        = "ng",
    concentration = "ng/mL"
  )

  covariateData <- list(
    DOSE = list(
      description        = paste(
        "Total administered acetaminophen dose per subject. Used by the",
        "enterohepatic-recirculation term F_abs * DOSE in the plasma ODE,",
        "which adds a continuous mass input to plasma proportional to the",
        "originally administered dose."
      ),
      units              = "ng",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "For the Westerhout 2012 rat protocol (15 mg/kg IV infusion over",
        "10 min in 225-275 g rats), DOSE = 15 mg/kg * WT_kg * 1e6 ng/mg.",
        "A 250-g rat receives DOSE = 3.75e6 ng. The same value must also",
        "be supplied via the rxode2 dose event amt column on the infusion",
        "row; DOSE here is a per-subject covariate carrying the originally",
        "administered total mass so the F_abs continuous input term can be",
        "computed inside model() without reading the event-table amt."
      ),
      source_name        = "DOSE"
    )
  )

  population <- list(
    species        = "rat (male Wistar WU)",
    n_subjects     = 24L,
    n_studies      = 1L,
    age_range      = "adult",
    weight_range   = "225-275 g (mean ~250 g)",
    sex_female_pct = 0,
    disease_state  = paste(
      "Healthy adult male Wistar WU rats (Charles River, Maastricht, NL).",
      "Housed under standard conditions (21 C, 60 percent humidity, 12/12 h",
      "light/dark, ad libitum food and water). 24 rats randomized to three",
      "two-probe microdialysis combinations (ST + LV, ST + CM, or LV + CM;",
      "n = 8 per combination) where ST = striatum (brain ECF), LV = lateral",
      "ventricle CSF, and CM = cisterna magna CSF. Cannulas implanted in",
      "left femoral artery (blood sampling) and vein (drug administration)",
      "under isoflurane anesthesia; CMA/12 microdialysis guides chronically",
      "implanted in the named brain regions. The study protocol was",
      "approved by the Animal Ethics Committee of Leiden University",
      "(UDEC no. 07068)."
    ),
    dose_range     = paste(
      "Single 15 mg/kg acetaminophen IV infusion over 10 min delivered via",
      "automated pump (Pump 22 Multiple Syringe Pump, Harvard Apparatus)",
      "at 200 uL/min/kg. A second 15 mg/kg infusion at t = 240 min was",
      "given solely to determine plasma protein binding at C_max; the",
      "structural PK model fits only the first-dose 0-240 min interval."
    ),
    sampling       = paste(
      "Blood samples at t = -5 (blank), 2, 7, 10, 15, 30, 60, 120, 180,",
      "and 240 min (100 uL each, from the arterial cannula). Microdialysate",
      "collected at 10-min intervals from t = -1 h to t = +2 h and at 20-min",
      "intervals from t = +2 h to t = +4 h. Acetaminophen quantified by",
      "HPLC with electrochemical detection."
    ),
    regions        = "Leiden University, The Netherlands",
    plasma_protein_binding = paste(
      "Linear plasma protein binding determined by Centrifree",
      "ultrafiltration: fu_p = 0.805 +/- 0.042 (i.e., 19.5 +/- 4.2 percent",
      "bound). The model was fitted to unbound plasma concentrations",
      "(C_pl_u = fu_p * C_pl_total) after data preprocessing; fu_p does",
      "not enter the structural ODEs as written in the Appendix of",
      "Westerhout 2012, so it is documented here as population metadata",
      "rather than as an ini() parameter."
    ),
    in_vivo_recovery = paste(
      "Microdialysis in-vivo recoveries determined by reverse dialysis in",
      "a separate cohort of 12 rats: striatum 12.0 +/- 3.3 percent, lateral",
      "ventricle 8.1 +/- 3.8 percent, cisterna magna 8.6 +/- 4.7 percent.",
      "Brain ECF and CSF concentrations reported in Fig. 1 of Westerhout",
      "2012 are dialysate concentrations divided by the location-specific",
      "in-vivo recovery; this conversion is upstream of the model."
    ),
    notes          = paste(
      "Recirculation: the additional plateau-shaping zero-order input",
      "F_abs * DOSE represents reabsorption of biliary acetaminophen-",
      "glucuronide and acetaminophen-sulfate metabolites that are",
      "hydrolyzed back to acetaminophen in the small intestine. Estimation",
      "method: NONMEM 6.2 (Icon Development Solutions); ruggedness checked",
      "by n = 100 nonparametric case-resampling bootstraps and visual",
      "predictive check (1000 simulated trials, 90 percent prediction",
      "interval)."
    )
  )

  ini({
    # =========================================================================
    # Westerhout 2012 Table I, rat column. Reported estimate +/- standard
    # error; SE in the trailing comment. Clearances converted from the
    # paper's heterogeneous units (CL10 / Q12 in mL/min; CL13 / CL31 /
    # CL14 / CL41 / CL16 / CL61 in uL/min) to a single mL/min basis for
    # internal unit consistency with the volumes in mL.
    # =========================================================================

    # Plasma elimination clearance, CL10 (Table I).
    lcl10 <- log(13.8)            ; label("Plasma elimination clearance CL10 (mL/min)")             # Table I, rat: 13.8 +/- 1.0 mL/min

    # Intercompartmental clearance to peripheral tissue, Q12 (Table I).
    lq12 <- log(45.1)             ; label("Plasma-peripheral inter-compartmental clearance Q12 (mL/min)")  # Table I, rat: 45.1 +/- 5.8 mL/min

    # Plasma to brain ECF clearance CL13 (paper reports 165 +/- 39 uL/min;
    # converted to mL/min = 165e-3 = 0.165).
    lcl13 <- log(0.165)           ; label("Plasma-to-brain-ECF clearance CL13 (mL/min)")            # Table I, rat: 165 +/- 39 uL/min = 0.165 mL/min

    # Brain ECF to plasma clearance CL31 (paper reports 198 +/- 24 uL/min).
    lcl31 <- log(0.198)           ; label("Brain-ECF-to-plasma clearance CL31 (mL/min)")            # Table I, rat: 198 +/- 24 uL/min = 0.198 mL/min

    # Plasma to CSF lateral-ventricle clearance CL14 (paper 2.9 +/- 1.3 uL/min).
    lcl14 <- log(2.9e-3)          ; label("Plasma-to-CSF-lateral-ventricle clearance CL14 (mL/min)")    # Table I, rat: 2.9 +/- 1.3 uL/min = 0.0029 mL/min

    # CSF lateral-ventricle to plasma clearance CL41 (paper 5.0 +/- 2.1 uL/min).
    lcl41 <- log(5.0e-3)          ; label("CSF-lateral-ventricle-to-plasma clearance CL41 (mL/min)")    # Table I, rat: 5.0 +/- 2.1 uL/min = 0.005 mL/min

    # Plasma to CSF cisterna-magna clearance CL16 (paper 0.8 +/- 0.4 uL/min).
    lcl16 <- log(0.8e-3)          ; label("Plasma-to-CSF-cisterna-magna clearance CL16 (mL/min)")       # Table I, rat: 0.8 +/- 0.4 uL/min = 0.0008 mL/min

    # CSF cisterna-magna to plasma clearance CL61 (paper 4.5 +/- 0.9 uL/min).
    lcl61 <- log(4.5e-3)          ; label("CSF-cisterna-magna-to-plasma clearance CL61 (mL/min)")       # Table I, rat: 4.5 +/- 0.9 uL/min = 0.0045 mL/min

    # Peripheral volume of distribution V_per (paper 188 +/- 11 mL).
    lvp <- log(188)               ; label("Peripheral volume of distribution V_per (mL)")           # Table I, rat: 188 +/- 11 mL

    # Enterohepatic recirculation rate F_abs (paper 0.025 percent / min = 2.5e-4 / min).
    # SE not reported in Table I for F_abs; paper text describes this as
    # estimated by the model fit ("This additional, infusion-like, drug
    # input represents the fraction F_abs ... that is recirculated within
    # the study duration (240 min)"). Encoded as estimated without
    # fixed() because the source describes it as a fit parameter.
    lfabs <- log(2.5e-4)          ; label("Enterohepatic recirculation rate constant F_abs (1/min)")  # Table I, rat: 0.025 percent/min = 2.5e-4 /min (SE not reported)

    # =========================================================================
    # IIV variances. Westerhout 2012 reports per-parameter eta values as
    # the NONMEM OMEGA diagonal (variance scale). Interpreted as
    # omega^2 such that the lognormal SD on the back-transformed scale is
    # sqrt(omega^2). Off-diagonal correlations are not reported in the
    # paper, so each eta enters as an uncorrelated single-parameter IIV.
    # IIVs are reported for CL10, CL13, CL14, and CL16 only; the
    # remaining parameters were fitted without IIV (set to zero).
    # =========================================================================
    etalcl10 ~ 0.03               # Table I, rat: eta_CL10 = 0.03 +/- 0.01 (variance)
    etalcl13 ~ 0.45               # Table I, rat: eta_CL13 = 0.45 +/- 0.25 (variance)
    etalcl14 ~ 0.28               # Table I, rat: eta_CL14 = 0.28 +/- 0.13 (variance)
    etalcl16 ~ 1.11               # Table I, rat: eta_CL16 = 1.11 +/- 0.54 (variance)

    # =========================================================================
    # Proportional residual error. Westerhout 2012 reports per-output
    # epsilon values; interpreted as the NONMEM SIGMA diagonal (variance
    # scale). Back-transformed to a proportional residual SD via
    # sqrt(sigma^2). Per the canonical multi-output residual-error
    # convention (parameter-names.md), the parent output Cc uses the
    # bare propSd and the named outputs use propSd_<output>.
    # =========================================================================
    propSd          <- sqrt(0.08); label("Proportional residual SD on plasma Cc (fraction)")        # Table I, rat: epsilon_plasma = 0.08 +/- 0.02 (variance); SD = sqrt(0.08)
    propSd_Cbrain_ecf <- sqrt(0.14); label("Proportional residual SD on brain ECF Cbrain_ecf (fraction)") # Table I, rat: epsilon_brain_ECF = 0.14 +/- 0.03 (variance); SD = sqrt(0.14)
    propSd_Ccsf_lv  <- sqrt(0.19); label("Proportional residual SD on CSF lateral-ventricle Ccsf_lv (fraction)") # Table I, rat: epsilon_CSF_LV = 0.19 +/- 0.05 (variance); SD = sqrt(0.19)
    propSd_Ccsf_cm  <- sqrt(0.18); label("Proportional residual SD on CSF cisterna-magna Ccsf_cm (fraction)") # Table I, rat: epsilon_CSF_CM = 0.18 +/- 0.04 (variance); SD = sqrt(0.18)
  })

  model({
    # =========================================================================
    # Rat physiological volumes (mL) and bulk fluid flow rates (mL/min)
    # fixed at literature values per Table I of Westerhout 2012. uL
    # values converted to mL by dividing by 1000.
    # =========================================================================
    V_pl   <- 10.6                              # Table I, rat: V_pl = 10.6 mL (ref 49 Lee 1985)
    V_ICS  <- 1440 / 1000                       # Table I, rat: V_ICS = 1440 uL = 1.44 mL (ref 23 Cserr 1965; 80 percent of 1.8 g brain)
    V_ECF_phys <- 290 / 1000                    # Table I, rat: V_ECF = 290 uL = 0.290 mL (ref 30 Abbott 2004)
    # Effective brain compartment volume: paper text page 5 states
    # "the brain ICS volume is added to the brain ECF volume to account
    # for the total brain volume" (because brain ICS concentrations are
    # assumed equal to brain ECF concentrations). V_brain therefore
    # equals V_ECF_phys + V_ICS = 0.290 + 1.44 = 1.73 mL.
    V_brain <- V_ECF_phys + V_ICS               # 1.73 mL effective brain distribution volume
    V_LV   <- 50  / 1000                        # Table I, rat: V_LV = 50 uL = 0.050 mL (refs 25-26 Condon 1986, Kohn 1991; 17 percent of 300 uL total CSF)
    V_TFV  <- 50  / 1000                        # Table I, rat: V_TFV = 50 uL = 0.050 mL (ref 29 Adam 1978; 17 percent unaccounted-for CSF, in 3rd + 4th ventricles)
    V_CM   <- 17  / 1000                        # Table I, rat: V_CM = 17 uL = 0.017 mL (refs 27-28 Robertson 1949, Adam 1978; 6 percent of total CSF)
    V_SAS  <- 180 / 1000                        # Table I, rat: V_SAS = 180 uL = 0.180 mL (refs 24, 29 Bass 1973, Adam 1978; 60 percent of total CSF)
    Q_ECF  <- 0.2  / 1000                       # Table I, rat: Q_ECF = 0.2 uL/min = 2e-4 mL/min (refs 30-31 Abbott 2004, Cserr 1981)
    Q_CSF  <- 2.2  / 1000                       # Table I, rat: Q_CSF = 2.2 uL/min = 2.2e-3 mL/min (ref 13 Cserr 1965)

    # =========================================================================
    # Individual structural parameters: back-transform from log scale,
    # apply IIV where the paper reports it.
    # =========================================================================
    cl10 <- exp(lcl10 + etalcl10)
    q12  <- exp(lq12)
    cl13 <- exp(lcl13 + etalcl13)
    cl31 <- exp(lcl31)
    cl14 <- exp(lcl14 + etalcl14)
    cl41 <- exp(lcl41)
    # Per paper Results: "Since we have no measurements of the
    # concentrations in the third and fourth ventricle, the transfer
    # clearance between plasma and third and fourth ventricle was assumed
    # to be equal to the transfer clearance between plasma and lateral
    # ventricle." CL15 and CL51 are therefore structurally tied to CL14
    # and CL41 (no separate parameter, IIV propagates implicitly via the
    # CL14 / CL41 etas).
    cl15 <- cl14
    cl51 <- cl41
    cl16 <- exp(lcl16 + etalcl16)
    cl61 <- exp(lcl61)
    vp   <- exp(lvp)
    fabs <- exp(lfabs)

    # =========================================================================
    # Micro-constants (Westerhout 2012 Appendix, "Where:" block on
    # page 551). Each k_ij = CL_ij / V_origin where the origin volume is
    # the compartment that the clearance is flowing OUT OF.
    # =========================================================================
    k10 <- cl10 / V_pl     # Appendix: k10 = CL10 / V_pl
    k12 <- q12  / V_pl     # Appendix: k12 = Q12  / V_pl
    k21 <- q12  / vp       # Appendix: k21 = Q12  / V_per
    k13 <- cl13 / V_pl     # Appendix: k13 = CL13 / V_pl
    k31 <- cl31 / V_brain  # Appendix: k31 = CL31 / V_ECF (text: V_ECF + V_ICS)
    k14 <- cl14 / V_pl     # Appendix: k14 = CL14 / V_pl
    k41 <- cl41 / V_LV     # Appendix: k41 = CL41 / V_LV
    k15 <- cl15 / V_pl     # Appendix: k15 = CL15 / V_pl
    k51 <- cl51 / V_TFV    # Appendix: k51 = CL51 / V_TFV
    k16 <- cl16 / V_pl     # Appendix: k16 = CL16 / V_pl
    k61 <- cl61 / V_CM     # Appendix: k61 = CL61 / V_CM

    # =========================================================================
    # Plasma compartment (Appendix "Plasma:" block). dose is the
    # event-driven IV infusion input; F_abs * DOSE is the continuous
    # enterohepatic-recirculation input proportional to the originally
    # administered total dose (paper Results: "additional, infusion-like
    # drug input ... fraction (F_abs) of the administered acetaminophen
    # that is recirculated within the study duration (240 min)").
    # =========================================================================
    d/dt(central) <-
      - k12 * central + k21 * peripheral1
      - k13 * central + k31 * brain_ecf
      - k14 * central + k41 * csf_lv
      - k15 * central + k51 * csf_tfv
      - k16 * central + k61 * csf_cm
      + (Q_CSF / V_SAS) * csf_sas
      - k10 * central
      + fabs * DOSE

    # Periphery (Appendix "Periphery:" block).
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Brain ECF (Appendix "Brain ECF:" block). The (Q_ECF / V_brain) * A
    # term carries the ECF-to-CSF bulk flow into CSF_LV (next equation).
    d/dt(brain_ecf) <- k13 * central - k31 * brain_ecf - (Q_ECF / V_brain) * brain_ecf

    # CSF lateral ventricle (Appendix "CSF_LV:" block). Receives the
    # brain-ECF bulk-flow input (Q_ECF / V_brain * A_ECF) and exports
    # downstream to CSF_TFV via Q_CSF.
    d/dt(csf_lv) <- k14 * central - k41 * csf_lv + (Q_ECF / V_brain) * brain_ecf - (Q_CSF / V_LV) * csf_lv

    # CSF third + fourth ventricle (Appendix "CSF_TFV:" block). Receives
    # Q_CSF * C_LV from CSF_LV and exports to CSF_CM via Q_CSF.
    d/dt(csf_tfv) <- k15 * central - k51 * csf_tfv + (Q_CSF / V_LV) * csf_lv - (Q_CSF / V_TFV) * csf_tfv

    # CSF cisterna magna (Appendix "CSF_CM:" block). Receives Q_CSF * C_TFV
    # from CSF_TFV and exports to CSF_SAS via Q_CSF.
    d/dt(csf_cm) <- k16 * central - k61 * csf_cm + (Q_CSF / V_TFV) * csf_tfv - (Q_CSF / V_CM) * csf_cm

    # CSF subarachnoid space (Appendix "CSF_SAS:" block). Receives
    # Q_CSF * C_CM from CSF_CM and exports back to plasma via Q_CSF (the
    # plasma compartment's "+ (Q_CSF / V_SAS) * A_SAS" return term above).
    d/dt(csf_sas) <- (Q_CSF / V_CM) * csf_cm - (Q_CSF / V_SAS) * csf_sas

    # =========================================================================
    # Observation variables and residual-error mapping. The fitted
    # outputs are unbound plasma (Cc), brain ECF (Cbrain_ecf), CSF
    # lateral-ventricle (Ccsf_lv), and CSF cisterna-magna (Ccsf_cm).
    # CSF third-fourth-ventricle (Ccsf_tfv) and subarachnoid (Ccsf_sas)
    # concentrations are predicted by the structural model but were not
    # measured in Westerhout 2012; they are exposed as derived outputs
    # without residual error so downstream users can plot them.
    # =========================================================================
    Cc          <- central     / V_pl
    Cbrain_ecf  <- brain_ecf   / V_brain
    Ccsf_lv     <- csf_lv      / V_LV
    Ccsf_tfv    <- csf_tfv     / V_TFV
    Ccsf_cm     <- csf_cm      / V_CM
    Ccsf_sas    <- csf_sas     / V_SAS

    Cc         ~ prop(propSd)
    Cbrain_ecf ~ prop(propSd_Cbrain_ecf)
    Ccsf_lv    ~ prop(propSd_Ccsf_lv)
    Ccsf_cm    ~ prop(propSd_Ccsf_cm)
  })
}
