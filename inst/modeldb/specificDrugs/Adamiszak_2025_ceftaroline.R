Adamiszak_2025_ceftaroline <- function() {
  description <- paste(
    "Five-compartment population PK model for ceftaroline in critically ill",
    "adults with and without continuous renal replacement therapy (CRRT): a",
    "two-compartment patient model (first-order elimination from the central,",
    "pre-filter compartment) joined to three CRRT-machine compartments (filter,",
    "post-filter return limb and effluent) whose transfer rate constants are",
    "circuit flow divided by compartment volume, all gated on by a time-varying",
    "CRRT-running indicator. Creatinine clearance (power model centred at 159",
    "mL/min/1.73 m^2) on clearance; CRRT modality (CVVHD, CVVHDF) on the",
    "central volume; any CRRT on the peripheral volume and intercompartmental",
    "clearance. Between-occasion variability on clearance, between-subject",
    "variability on the central and peripheral volumes and intercompartmental",
    "clearance, and proportional residual error for each of the three sampled",
    "matrices (systemic pre-filter plasma, post-filter plasma, effluent)."
  )
  reference <- paste(
    "Adamiszak A, Pietrzkiewicz K, Bartkowska-Sniatkowska A, Smuszkiewicz P,",
    "Kusza K, Grzeskowiak E, Bienert A (2025).",
    "Do Critically Ill Patients Undergoing Continuous Renal Replacement Therapy",
    "Require Ceftaroline Dosage Adjustments? Ceftaroline PopPK Model and Dosage",
    "Simulations with the Probability of Target Attainment Analysis Based on",
    "Retrospective Data.",
    "Antibiotics 14(4):347.",
    "doi:10.3390/antibiotics14040347."
  )
  vignette <- "Adamiszak_2025_ceftaroline"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # The structure is the supplement's Monolix (mlxtran) final-model syntax,
  # which is identical to ModelFile/Syntax_v2_SameClbody+IS_CRRT+Q.txt in the
  # authors' Simulx project (github.com/arkadiusz-adamiszak/
  # Ceftaroline_ICU_Retrospective). Paper -> canonical state names:
  #   Acent  -> central            (Figure 1 'Central Pre-filter')
  #   Aper   -> peripheral1
  #   Acrrt  -> circuit_plasma     (Figure 1 'CRRT'; blood side of the filter,
  #                                 red-cell phase lumped in, as the paper does)
  #   Apost  -> circuit_postfilter (Figure 1 'Post-filter')
  #   Aeff   -> circuit_dialysate  (Figure 1 'Effluent'; a flow-through volume
  #                                 that gains drug from the filter and drains
  #                                 onward at Qeff, i.e. the filtrate/dialysate
  #                                 side, NOT the accumulating effluent bag
  #                                 that circuit_effluent denotes)
  compartmentData <- list(
    central = list(analyte = "ceftaroline", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "ceftaroline", units = "mg", specimen = "plasma", verified = TRUE),
    circuit_plasma = list(analyte = "ceftaroline", units = "mg", specimen = "whole blood", verified = TRUE),
    circuit_postfilter = list(analyte = "ceftaroline", units = "mg", specimen = "plasma", verified = TRUE),
    circuit_dialysate = list(analyte = "ceftaroline", units = "mg", specimen = "dialysate", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description = "Creatinine clearance by the MDRD formula, BSA-normalised; power-model covariate on clearance centred at the analysis-population median 159 mL/min/1.73 m^2",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Table 1: median 159 [range 42-309] mL/min/1.73 m^2, estimated by the",
        "Modification of Diet in Renal Disease (MDRD) formula. Methods 4.2:",
        "continuous covariates enter as theta_pop * (COV / median COV)^theta_cov;",
        "the centring value 159 is the Table 1 median and is confirmed by the",
        "Simulx project's covariate transform logtClCr_Pros = log(ClCr_Pros/159).",
        "Occasion-varying: for the 18 patients of one source study two CrCL",
        "values were available and were carried per occasion together with the",
        "between-occasion variability on clearance (Methods 4.2). For the CRRT",
        "simulations the paper assumes CrCL 58 mL/min/1.73 m^2, the median of",
        "the four CRRT patients (Methods 4.3)."
      ),
      source_name = "ClCr_Pros"
    ),
    RRT_CVVHD_STATUS = list(
      description = "Continuous venovenous hemodialysis (CVVHD) modality indicator; log-additive effect on the central volume",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (not on CVVHD); with RRT_CVVHDF_STATUS also 0 the subject is the no-CRRT reference",
      notes = paste(
        "Decomposed from the three-level categorical CRRT modality covariate",
        "(no CRRT / CVVHD / CVVHDF, reference no CRRT; Simulx covariate",
        "tCRRT_NaRef_CVVHD_CVVHDF built on the data column 'CRRT type'). Subject",
        "level and time-fixed: modality is a Monolix covariate on the individual",
        "parameters, separate from the time-varying CRRT-running regressor",
        "IS_CRRT (RRT_CRRT_ACTIVE). 2 of the 4 CRRT patients (Table 1). The sum",
        "RRT_CVVHD_STATUS + RRT_CVVHDF_STATUS is the pooled any-CRRT indicator",
        "that the paper places on the peripheral volume and Q, so",
        "RRT_CRRT_STATUS is not carried as a separate column."
      ),
      source_name = "CRRT type"
    ),
    RRT_CVVHDF_STATUS = list(
      description = "Continuous venovenous hemodiafiltration (CVVHDF) modality indicator; log-additive effect on the central volume",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (not on CVVHDF); with RRT_CVVHD_STATUS also 0 the subject is the no-CRRT reference",
      notes = paste(
        "Second non-reference level of the three-level CRRT modality covariate;",
        "the reference is NO CRRT, not CVVHD (contrast Broeker_2018_tigecycline.R,",
        "whose cohort is all-CRRT and uses CVVHD as reference). Subject level and",
        "time-fixed. 2 of the 4 CRRT patients (Table 1)."
      ),
      source_name = "CRRT type"
    ),
    RRT_CRRT_ACTIVE = list(
      description = "Time-varying CRRT-running indicator; switches every circuit flow on (1) or off (0)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (circuit not running)",
      notes = paste(
        "The mlxtran regressor IS_CRRT, which multiplies every circuit term of",
        "the ODE system. Setting it to 0 freezes the three machine compartments",
        "and reduces the model to a two-compartment patient model; the paper's",
        "Figure 5 switches it off for 24 h mid-course and back on with a",
        "different effluent flow rate (Simulx regressor 'Real_Time', times",
        "0 / 48 / 72 h). It does NOT switch the modality effects on the",
        "volumes and Q, which are subject-level."
      ),
      source_name = "IS_CRRT"
    ),
    BFR = list(
      description = "CRRT blood flow rate; drives filter inflow and post-filter return as Qb / V",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Table 1: median 250 [range 250-300] mL/min. The paper's model and",
        "dataset carry Qb in L/h (15 L/h = 250 mL/min, Methods 4.3); the",
        "conversion x 60 / 1000 is applied inside model(). Meaningful only",
        "while RRT_CRRT_ACTIVE = 1; set to 0 for patients not on CRRT (Figure 1",
        "footnote: 'When a patient is not under CRRT, all flow rates are set to",
        "zero')."
      ),
      source_name = "Qb"
    ),
    RRT_CRRT_EFFLUENT_FLOW = list(
      description = "CRRT total effluent flow rate; drives filter-to-effluent transfer and effluent washout as Qeff / V",
      units = "mL/h",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Table 1: median 3200 [range 2700-3660] mL/h. The paper's model and",
        "dataset carry Qeff in L/h (simulations at 2.5, 3.0, 3.5 and 4.0 L/h,",
        "Methods 4.3); the conversion / 1000 is applied inside model().",
        "Meaningful only while RRT_CRRT_ACTIVE = 1; set to 0 otherwise."
      ),
      source_name = "Qeff"
    ),
    OCC = list(
      description = "Occasion index for the between-occasion variability on clearance",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = paste(
        "Up to two occasions per patient: 15 patients of one source study",
        "contribute two occasions in the authors' deposited dataset, each with",
        "its own CrCL (Methods 4.2 states 18 patients had two CrCL",
        "measurements). Decomposed inside model() into the indicators oc1 / oc2",
        "selecting the per-occasion eta, because rxode2 parses but cannot",
        "simulate the `eta ~ var | occ` syntax; the occasion-2 variance is",
        "fixed() to occasion 1 (one IOV magnitude reported). Records with any",
        "other OCC value receive no IOV."
      ),
      source_name = "OCC"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 33L,
    n_studies = 4L,
    age_range = "21-77 years (median 51)",
    age_median = "51 years",
    weight_range = "49-111.8 kg (median 75.5)",
    weight_median = "75.5 kg",
    sex_female_pct = 36,
    race_ethnicity = "Not reported",
    disease_state = paste(
      "Critically ill adults in intensive care treated with ceftaroline",
      "fosamil; 4 of 33 (12%) on continuous renal replacement therapy",
      "(2 CVVHD, 2 CVVHDF; blood flow 250 [250-300] mL/min, effluent flow",
      "3200 [2700-3660] mL/h) (Table 1)."
    ),
    renal_function = "MDRD creatinine clearance median 159 [42-309] mL/min/1.73 m^2 (Table 1); the four CRRT patients 42-67 mL/min/1.73 m^2.",
    dose_range = "400-600 mg ceftaroline fosamil every 8-12 h as a 1 h intravenous infusion (Methods 4.1).",
    regions = "Pooled from four published studies (Alarcia-Lacalle 2021, Chauzy 2019, Chauzy 2022, Kalaria 2021)",
    notes = paste(
      "Concentrations were DIGITISED (WebPlotDigitizer 4.7) from figures of the",
      "four source publications rather than taken from patient-level records:",
      "296 plasma concentrations from 29 non-CRRT patients, plus 24 pre-filter",
      "(systemic), 23 post-filter and 23 effluent concentrations from the 4",
      "CRRT patients (Results 2.1). Eight BLQ values handled with the M3",
      "method. Fit with SAEM in Monolix 2024R1."
    )
  )

  ini({
    # Structural parameters -- Table 2, 'Final Model' column (values recovered
    # from the PDF layout; the article's HTML table is mis-rendered).
    lcl <- log(11.23)
    label("Clearance from the central compartment at CrCL 159 mL/min/1.73 m^2, Cl (L/h)") # Table 2 'Cl (L/h)' 11.23 (RSE 5.7%)
    e_crcl_cl <- 0.68
    label("Power exponent of CrCL/159 on clearance (unitless)") # Table 2 'ClCr on Cl' 0.68 (RSE 16.1%)

    lvc <- log(20.04)
    label("Central (pre-filter) volume without CRRT, VCentral (L)") # Table 2 'VCentral (L)' 20.04 (RSE 4.1%)
    e_cvvhd_vc <- 0.47
    label("Log-scale CVVHD effect on the central volume (unitless)") # Table 2 'CVVHD on VCentral' 0.47 (RSE 18.2%)
    e_cvvhdf_vc <- -0.42
    label("Log-scale CVVHDF effect on the central volume (unitless)") # Table 2 'CVVHDF on VCentral' -0.42 (RSE 21.3%)

    lq <- log(8.80)
    label("Intercompartmental clearance without CRRT, Q (L/h)") # Table 2 'Q (L/h)' 8.80 (RSE 20.7%)
    e_crrt_q <- 0.99
    label("Log-scale any-CRRT effect on Q (unitless)") # Table 2 'CRRT on Q' 0.99 (RSE 26.6%)

    lvp <- log(16.85)
    label("Peripheral volume without CRRT, VPeripheral (L)") # Table 2 'VPeripheral (L)' 16.85 (RSE 10.5%)
    e_crrt_vp <- 0.99
    label("Log-scale any-CRRT effect on the peripheral volume (unitless)") # Table 2 'CRRT on VPeripheral' 0.99 (RSE 30.5%)

    lvfilter <- log(2.69)
    label("CRRT filter (blood-side) compartment volume, VCRRT (L)") # Table 2 'VCRRT (L)' 2.69 (RSE 42.3%)
    lvpostfilter <- log(1.15)
    label("CRRT post-filter compartment volume, VPost-filter (L)") # Table 2 'VPost-filter (L)' 1.15 (RSE 3.5%)
    lvdialysate <- log(0.60)
    label("CRRT effluent compartment volume, VEffluent (L)") # Table 2 'VEffluent (L)' 0.60 (RSE 6.1%)

    # Random effects -- Table 2 reports each as %CV; converted to log-normal
    # variance omega^2 = log(1 + CV^2). No IIV on clearance: only IOV
    # (Results 2.2, 'Due to the instability of the model').
    etaiov_cl_1 ~ 0.115309 # Table 2 'IOV on Cl (%CV)' 34.96 (RSE 12.3%) -> log(1 + 0.3496^2)
    etaiov_cl_2 ~ fixed(0.115309) # shared magnitude with occasion 1 (one IOV value reported)
    etalvc ~ 0.00136806 # Table 2 'IIV on VCentral (%CV)' 3.70 (RSE 57.8%) -> log(1 + 0.0370^2)
    etalq ~ 0.0265407 # Table 2 'IIV on Q (%CV)' 16.40 (RSE 58.5%) -> log(1 + 0.1640^2)
    etalvp ~ 0.138889 # Table 2 'IIV on VPeripheral (%CV)' 38.60 (RSE 20.8%) -> log(1 + 0.3860^2)

    # Residual error -- proportional for each matrix (Results 2.2; mlxtran
    # errorModel=proportional).
    propSd <- 0.29
    label("Proportional residual SD, systemic (pre-filter) plasma, bPre (fraction)") # Table 2 'bPre' 0.29 (RSE 5.0%)
    propSd_Cpostfilter <- 0.16
    label("Proportional residual SD, post-filter plasma, bPost (fraction)") # Table 2 'bPost' 0.16 (RSE 21.0%)
    propSd_Ceffluent <- 0.13
    label("Proportional residual SD, effluent, bEff (fraction)") # Table 2 'bEff' 0.13 (RSE 16.9%)
  })
  model({
    # Occasion indicators for the between-occasion variability on clearance.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    iov_cl <- oc1 * etaiov_cl_1 + oc2 * etaiov_cl_2

    # Pooled any-CRRT indicator (Simulx tCRRT_NaRef_CRRT: CVVHD and CVVHDF in
    # one group against no CRRT).
    crrt <- RRT_CVVHD_STATUS + RRT_CVVHDF_STATUS

    # Individual parameters (Methods 4.2: continuous covariates as
    # (COV / median)^theta, categorical as exp(theta)).
    cl <- exp(lcl + iov_cl) * (CRCL / 159)^e_crcl_cl
    vc <- exp(lvc + e_cvvhd_vc * RRT_CVVHD_STATUS + e_cvvhdf_vc * RRT_CVVHDF_STATUS + etalvc)
    q <- exp(lq + e_crrt_q * crrt + etalq)
    vp <- exp(lvp + e_crrt_vp * crrt + etalvp)
    vfilter <- exp(lvfilter)
    vpostfilter <- exp(lvpostfilter)
    vdialysate <- exp(lvdialysate)

    # Circuit flows in L/h: BFR is charted in mL/min, the effluent flow in mL/h.
    qb <- BFR * 60 / 1000
    qeff <- RRT_CRRT_EFFLUENT_FLOW / 1000

    # Rate constants, k = flow / volume (Methods 4.2 transfer formula and the
    # mlxtran EQUATION block).
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    k_central_filter <- qb / vc # mlxtran k13
    k_postfilter_central <- qb / vpostfilter # mlxtran k41
    k_filter_postfilter <- qb / vfilter # mlxtran k34
    k_filter_dialysate <- qeff / vfilter # mlxtran k35
    k_dialysate_out <- qeff / vdialysate # mlxtran kelcrrt

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1 +
      RRT_CRRT_ACTIVE * (-k_central_filter * central + k_postfilter_central * circuit_postfilter)
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(circuit_plasma) <- RRT_CRRT_ACTIVE *
      (k_central_filter * central - k_filter_dialysate * circuit_plasma - k_filter_postfilter * circuit_plasma)
    d/dt(circuit_postfilter) <- RRT_CRRT_ACTIVE *
      (k_filter_postfilter * circuit_plasma - k_postfilter_central * circuit_postfilter)
    d/dt(circuit_dialysate) <- RRT_CRRT_ACTIVE *
      (k_filter_dialysate * circuit_plasma - k_dialysate_out * circuit_dialysate)

    Cc <- central / vc
    Cpostfilter <- circuit_postfilter / vpostfilter
    Ceffluent <- circuit_dialysate / vdialysate

    Cc ~ prop(propSd)
    Cpostfilter ~ prop(propSd_Cpostfilter)
    Ceffluent ~ prop(propSd_Ceffluent)
  })
}
