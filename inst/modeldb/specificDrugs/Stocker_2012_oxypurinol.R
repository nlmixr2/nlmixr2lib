Stocker_2012_oxypurinol <- function() {
  description <- "One-compartment population PK model for oxypurinol (the active metabolite of allopurinol) in adults with gout (Stocker 2012). First-order formation from allopurinol (Kfm taken as the apparent first-order absorption rate into the central compartment), one-compartment distribution, and first-order elimination. Apparent clearance (CL/Fm) is modified by raw Cockcroft-Gault creatinine clearance based on lean body weight (CRCL), concomitant any-class diuretic use (CONMED_DIUR; thiazide / furosemide / spironolactone pooled), and concomitant probenecid use (CONMED_PROBENECID), each via a linear-deviation multiplicative factor. Apparent volume (V/Fm) is allometrically scaled on lean body weight (LBW) with the volume exponent held fixed at the theoretical value of 1.0. The dose entered into the model is the oxypurinol-equivalent dose, taken as 0.9 x allopurinol dose per the paper's prior published assumption."
  reference   <- "Stocker SL, McLachlan AJ, Savic RM, Kirkpatrick CM, Graham GG, Williams KM, Day RO. The pharmacokinetics of oxypurinol in people with gout. Br J Clin Pharmacol. 2012 Sep;74(3):477-489. doi:10.1111/j.1365-2125.2012.04207.x"
  vignette    <- "Stocker_2012_oxypurinol"
  units       <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    CRCL = list(
      description        = "Creatinine clearance estimated by the Cockcroft-Gault equation using lean body weight as the size descriptor (NOT BSA-normalised).",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Stocker 2012 Methods (Covariate model development): 'creatinine clearance (CLCr) estimated using the Cockcroft-Gault equation corrected for LBW or using TBW'. The final model retained the LBW-based form. Enters CL/Fm as a linear-deviation multiplier centred on the cohort median 37.6 mL/min; not as a power scaling. Per-cohort reference value 37.6 mL/min reflects elderly gouty patients with substantial mild-to-moderate renal impairment (CLCr range 6.0-130.4 mL/min). Stored under the canonical CRCL register entry with the explicit non-BSA-normalised Cockcroft-Gault assay form documented here (matches Delattre_2010_amikacin.R precedent for raw mL/min Cockcroft-Gault).",
      source_name        = "CLCR"
    ),
    LBM = list(
      description        = "Lean body weight, paper notation LBW.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Stocker 2012 Methods cites reference [22] for the LBW formula (the James 1976 formula commonly used in the Sydney popPK group). Used as the size descriptor for V/Fm allometric scaling with the volume exponent held fixed at the theoretical value of 1.0 and a reference value of 60 kg (Stocker 2012 Table 3 footnote: 'For a typical patient who has a lean body weight of 60 kg'). LBW is also the size descriptor used to estimate CRCL via Cockcroft-Gault.",
      source_name        = "LBW"
    ),
    CONMED_DIUR = list(
      description        = "Concomitant any-class diuretic indicator (composite: thiazide + loop diuretic + potassium-sparing diuretic).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant diuretic of any class)",
      notes              = "Stocker 2012 Table 1 footnote ' Including furosemide, thiazide diuretics and spironolactone'. 46% of the cohort (72 of 155 gouty patients) were on a concomitant diuretic. Linear-deviation multiplicative effect on CL/Fm: cl *= (1 + (-0.294) * CONMED_DIUR), i.e. 29.4% lower apparent oxypurinol clearance when a diuretic is coadministered. Time-varying per occasion in the source dataset because changes in concomitant medication define a new occasion.",
      source_name        = "DIUR"
    ),
    CONMED_PROBENECID = list(
      description        = "Concomitant probenecid indicator.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant probenecid)",
      notes              = "Stocker 2012 Table 1: 13% of the cohort (20 of 155 gouty patients) were on chronic oral probenecid. Linear-deviation multiplicative effect on CL/Fm: cl *= (1 + 0.383 * CONMED_PROBENECID), i.e. 38.3% higher apparent oxypurinol clearance when probenecid is coadministered. Time-varying per occasion in the source dataset because changes in concomitant medication define a new occasion. The mechanism is probenecid inhibiting renal tubular reabsorption of oxypurinol (Stocker 2012 Discussion, citing references [21, 30, 38]).",
      source_name        = "PROB"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 155,
    n_studies       = 1,
    age_range       = "28.0 to 93.6 years",
    age_median      = "69.0 years",
    weight_range    = "42.5 to 139.0 kg total body weight",
    weight_median   = "83.0 kg total body weight",
    sex_female_pct  = 14.8,
    race_ethnicity  = "Not reported; cohort drawn from an Australian community and hospital setting in Sydney.",
    disease_state   = "Gout (n = 146) or asymptomatic hyperuricaemia (n = 9); all on chronic allopurinol for at least 7 days; 13% with tophi present; 28% (n = 44) with BMI > 30 kg/m^2.",
    dose_range      = "Allopurinol 50-400 mg/day oral (median 300 mg/day). Oxypurinol-equivalent doses entered into the model are 90% of the allopurinol dose per Stocker 2012 Methods (citing reference [6]).",
    regions         = "Australia (St Vincent's Hospital and community patients in New South Wales).",
    n_observations  = "1013 plasma oxypurinol concentrations (1-30 samples per patient, primarily across the 0-24 h dosing interval). Concentration range 0.14-70.9 mg/L. 23 samples below the limit of quantification (<2 mg/L) but above the limit of detection were included.",
    occasions       = "1-5 PK occasions per patient; occasion defined as a hospital readmission or a change in concomitant medication.",
    co_medication   = "All patients on at least one other medication (median 9, range 1-18). Significant covariates retained: concomitant diuretic (any class) and concomitant probenecid. Concomitant medications tested without retained effect: low-dose aspirin, colchicine, warfarin, angiotensin II receptor antagonists, beta-adrenoceptor antagonists, ACE inhibitors, calcium channel blockers, HMG-CoA reductase inhibitors, and antibiotics.",
    excluded        = "7 patients excluded prior to modelling: 4 not at steady-state, 2 undergoing dialysis, 1 in acute renal failure.",
    bov_cl_note     = "Source paper reports between-occasion variability (BOV) on CL/Fm of 21% CV (Stocker 2012 Table 3, 'BOV CL/Fm'). nlmixr2lib's ini() does not natively encode occasion-level random effects; the BSV CL/Fm reported in the model file (28% CV) does not include this BOV term. Users who need an occasion-aware simulation must add an additional occasion-level eta on lcl with variance log(1 + 0.21^2) = 0.04314 outside the model file.",
    notes           = "Cohort demographics in Stocker 2012 Table 1. Pharmacogenetic SNPs (SLC2A9, ABCG2, SLC22A13, SLC17A1, SLC22A11, SLC22A8) were tested but not retained after adjustment for CLCr, diuretics, and probenecid (Stocker 2012 Discussion). Bootstrap with 1000 replicates was used to refine parameter uncertainty (Stocker 2012 Table 3)."
  )

  ini({
    # Structural parameters - final-model NONMEM estimates (Stocker 2012
    # Table 3 'Covariate (final) model' column). All terms are apparent
    # because the study used oral allopurinol and no IV oxypurinol arm;
    # bioavailability is not separately identifiable from oxypurinol
    # formation, hence the paper's Fm notation.

    lka <- log(0.447); label("Apparent first-order absorption / formation rate constant Kfm (1/h)")        # Stocker 2012 Table 3: Kfm = 0.447; bootstrap median 0.438 (0.301, 0.628)
    lcl <- log(0.595); label("Apparent oral clearance CL/Fm at reference CRCL = 37.6 mL/min (L/h)")        # Stocker 2012 Table 3: CL/Fm = 0.595; bootstrap median 0.594 (0.548, 0.646); typical-patient CLCr 37.6 mL/min footnote
    lvc <- log(38.1);  label("Apparent central volume V/Fm at reference LBW = 60 kg (L)")                  # Stocker 2012 Table 3: V/Fm = 38.1; bootstrap median 38.3 (33.2, 44.4); typical-patient LBW 60 kg footnote

    # Allometric LBW exponent on V/Fm. Source paper writes the model
    # TVV/Fm = theta2 * (LBW/60)^theta9 (Stocker 2012 Results, equation
    # set on page 480) but does NOT report theta9 in Table 3 of the
    # final model. The paper also notes that the LBW scaling 'did not
    # lead to a significant improvement in the likelihood or decrease in
    # the BSV' and was retained for biological plausibility. The
    # canonical reading is that theta9 was held fixed at the theoretical
    # value 1.0 (linear scaling), consistent with the absence of a
    # Table 3 row for theta9 and with the typical fixed-exponent
    # convention for V/F allometric scaling (Anderson & Holford).
    e_lbm_vc <- fixed(1.0); label("Allometric LBW exponent on V/Fm (unitless, fixed)")                     # Stocker 2012 Results equation page 480: TVV/Fm = theta2 * (LBW/60)^theta9; theta9 not reported in Table 3, inferred fixed at theoretical 1.0

    # Covariate effects on CL/Fm. Linear-deviation multiplicative form
    # (NONMEM "additive on theta1" / linear-deviation scaling): CL =
    # TVCL * (1 + theta6 * (CRCL - 37.6)) * (1 + theta7 * DIUR) *
    # (1 + theta8 * PROB).
    e_crcl_cl              <-  0.0250; label("CRCL slope on CL/Fm: per mL/min of (CRCL - 37.6) deviation") # Stocker 2012 Table 3: theta6 = 0.0250 (0.021, 0.028)
    e_conmed_diur_cl       <- -0.294;  label("Fractional change in CL/Fm with concomitant any-class diuretic (CONMED_DIUR = 1)")    # Stocker 2012 Table 3: theta7 = -0.294 (-0.386, -0.207); -29.4%
    e_conmed_probenecid_cl <-  0.383;  label("Fractional change in CL/Fm with concomitant probenecid (CONMED_PROBENECID = 1)")      # Stocker 2012 Table 3: theta8 = 0.383 (0.264, 0.499); +38.3%

    # IIV. Stocker 2012 Methods (page 478) describes the IIV as
    # population parameter variability comprising between-subject
    # variability (BSV) plus between-occasion variability (BOV). Table 3
    # reports BSV CL/Fm = 28% CV, BSV V/Fm = 45% CV, and BOV CL/Fm =
    # 21% CV, with a positive correlation r = 0.33 between the CL and V
    # subject-level random effects. nlmixr2lib's ini() encodes the BSV
    # block here; the BOV is documented in population$bov_cl_note and
    # is the user's responsibility to add as a per-occasion eta if
    # needed.
    # CV% -> log-normal variance: omega^2 = log(CV^2 + 1)
    #   omega^2_cl = log(1 + 0.28^2) = 0.0755
    #   omega^2_vc = log(1 + 0.45^2) = 0.1844
    #   cov_cl_vc  = 0.33 * sqrt(0.0755 * 0.1844) = 0.0389
    etalcl + etalvc ~ c(0.0755, 0.0389, 0.1844)  # Stocker 2012 Table 3: BSV CL/Fm 28% CV, BSV V/Fm 45% CV, COV(CL, V) r = 0.33

    # Residual error. Stocker 2012 Methods describe a combined
    # additive + proportional residual error model retained as the
    # final-model RUV (Table 3 'Residual Error' block).
    addSd  <- 1.32;  label("Additive residual error (mg/L)")             # Stocker 2012 Table 3: additive 1.32; bootstrap median 1.32 (0.77, 1.64)
    propSd <- 0.131; label("Proportional residual error (fraction)")     # Stocker 2012 Table 3: proportional CV 13.1%; bootstrap median 12.9% (9.4%, 17.2%)
  })

  model({
    # Reference values used by the final covariate model (Stocker 2012
    # Table 3 footnotes).
    ref_crcl <- 37.6
    ref_lbm  <- 60

    # Linear-deviation covariate multiplier on CL/Fm. Combines the
    # three retained covariates by multiplication: the (1 + theta *
    # covariate) form preserves the meaning that each theta is the
    # fractional change in CL/Fm at the covariate's reference state
    # (CRCL = ref_crcl, DIUR = 0, PROB = 0).
    cov_cl <- (1 + e_crcl_cl              * (CRCL - ref_crcl)) *
              (1 + e_conmed_diur_cl       * CONMED_DIUR) *
              (1 + e_conmed_probenecid_cl * CONMED_PROBENECID)

    # Individual parameters.
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * cov_cl
    vc <- exp(lvc + etalvc) * (LBM / ref_lbm)^e_lbm_vc
    kel <- cl / vc

    # ODE system: one-compartment with first-order absorption /
    # oxypurinol formation. The depot compartment carries the
    # oxypurinol-equivalent dose (taken as 0.9 x allopurinol dose per
    # Stocker 2012 Methods) and the formation rate constant Kfm drives
    # appearance into the central compartment.
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Observation. Dose in mg, volume in L -> mg/L.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
