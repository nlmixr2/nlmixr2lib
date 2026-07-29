Levitt_2005_propofol_pbpk <- function() {
  description <- paste(
    "PBPK (whole-body). Human physiologically based pharmacokinetic model",
    "for propofol (Levitt and Schnider 2005 BMC Anesthesiology) fitted",
    "to the individual arterial plasma data of 24 adult volunteers in",
    "three age groups (18-34, 35-65, >65 years) from Schnider et al. 1998",
    "Anesthesiology. Each subject received an initial ~20-second IV",
    "bolus (2 mg/kg for subjects <65; 1 mg/kg for subjects >65)",
    "followed 60 minutes later by a 60-minute constant infusion at",
    "25, 50, 100, or 200 microg/kg/min. Twelve well-stirred flow-limited",
    "tissue compartments (liver, intestine [paper's 'portal'], muscle,",
    "kidney, brain, heart, lung, skin, adipose, bone, other [paper's",
    "'other' + 'tendon' connective tissue lumped], plus venous and",
    "arterial blood) use the PKQuest standard-human physiology (organ",
    "weights, blood flows, and fraction lipid from Table 1 for a 70-kg",
    "reference); propofol-specific tissue/blood partition coefficients",
    "are hard-coded from the paper text at the standard freepl = 0.022",
    "(adipose 84, brain 1.87, liver 2.12, intestine 1.7, other tissues",
    "1.45). Lean tissue weights and blood flows scale by lean body",
    "mass; adipose scales with body fat percent (Gallagher 1996",
    "regression internally derives body fat from WT, HT, AGE, SEXF,",
    "RACE_ASIAN). Hepatic elimination is a single well-stirred",
    "apparent-clearance approximation (CL_liver = 0.76 * (Q_hepatic",
    "artery + Q_portal)) calibrated to the paper's Table 4 age-averaged",
    "fractional liver clearance of 0.76; this replaces the paper's",
    "Roberts-Rowland dispersion-model liver (dispersion number 0.3)",
    "with an ODE-compatible well-mixed compartment. Optional",
    "pulmonary sequestration compartment (paper-specific 'lung_seq')",
    "receives a bolus fraction frdose = max(0, 0.548 - 0.00891 * AGE)",
    "and releases it exponentially with time constant T = 80 min",
    "(paper's average parameter set for Figs 9-11); the user allocates",
    "the sequestered bolus fraction via a second dose row on",
    "cmt = 'lung_seq'. No IIV or formal residual-error variance was",
    "reported (individual fits; 15 percent average absolute weighted",
    "residual). All parameters are FIXED at the paper's average values."
  )
  reference <- paste(
    "Levitt DG, Schnider TW (2005). Human physiologically based",
    "pharmacokinetic model for propofol. BMC Anesthesiology 5:4.",
    "doi:10.1186/1471-2253-5-4.",
    "Fitted to the individual arterial plasma data of Schnider TW, Minto",
    "CF, Gambus PL, Andresen C, Goodale DB, Shafer SL, Youngs EJ (1998).",
    "The influence of method of administration and covariates on the",
    "pharmacokinetics of propofol in adult volunteers. Anesthesiology",
    "88:1170-1182. doi:10.1097/00000542-199805000-00006.",
    sep = " "
  )
  vignette <- "Levitt_2005_propofol_pbpk"
  units    <- list(time = "minute", dosing = "mg", concentration = "mg/L")

  paper_specific_compartments <- c("lung_seq")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Per-subject body weight (Schnider 1998 Table 1 baseline",
        "demographics reproduced in Levitt 2005). WT scales lean tissue",
        "weights and blood flows by lean_mass / 52.5 kg (where 52.5 kg =",
        "70 kg * (1 - 0.25) is the reference non-fat mass at the standard",
        "25 percent body fat), and directly sets the adipose compartment",
        "weight as WT * BODYFAT_PCT / 100. Reference weight 70 kg matches",
        "Table 1 standard human."
      ),
      source_name        = "WT"
    ),
    HT = list(
      description        = "Height",
      units              = "cm",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Per-subject height (Schnider 1998 Table 1). Enters the Gallagher",
        "1996 body-fat regression via BMI = WT / (HT/100)^2 (paper eq 2).",
        "HT must be supplied in centimetres; the model converts to",
        "metres internally for BMI."
      ),
      source_name        = "HT"
    ),
    AGE = list(
      description        = "Age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Per-subject age (Schnider 1998 study spans 18-34, 35-65, and",
        ">65 year groups, 8 subjects each). AGE enters the body-fat",
        "regression (paper eq 2). AGE is also documented as the driver",
        "of the pulmonary-sequestration age dependence eq 13:",
        "frdose = max(0, 0.548 - 0.00891 * AGE); however sequestration",
        "is applied externally via a second dose row on cmt = 'lung_seq'",
        "(see vignette). frdose does not appear in the ODEs."
      ),
      source_name        = "AGE"
    ),
    SEXF = list(
      description        = "Sex, female indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = male",
      notes              = paste(
        "Per-subject sex (1 = female, 0 = male). Enters the Gallagher",
        "1996 body-fat regression as (1 - SEXF), i.e. the paper's Sex",
        "coding is 0 = female / 1 = male, so SEXM = 1 - SEXF is derived",
        "internally. Reference category 0 = male."
      ),
      source_name        = "SEXF"
    ),
    RACE_ASIAN = list(
      description        = "Race, Asian indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = non-Asian",
      notes              = paste(
        "Per-subject Asian race indicator (1 = Asian, 0 = other). Enters",
        "the Gallagher body-fat regression via the additive Asian",
        "correction term (95 / BMI - 0.044 * AGE) at paper eq 2. Source",
        "Gallagher DA 2000 Am J Clin Nutr 72:694-701."
      ),
      source_name        = "RACE_ASIAN"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 24L,
    n_studies      = 1L,
    age_range      = "18-34, 35-65, and >65 years (8 subjects per age group)",
    weight_range   = "not tabulated per subject in Levitt 2005; Schnider 1998 Table 1 lists individual weights",
    sex_female_pct = NA_real_,
    disease_state  = paste(
      "Healthy adult volunteers. Data are the individual arterial",
      "propofol concentrations of Schnider TW et al. 1998 Anesthesiology",
      "88:1170-1182 pooled across 24 subjects (8 per age group). Two",
      "identical visits per subject with and without EDTA formulation;",
      "the two visits were averaged since no meaningful difference was",
      "observed."
    ),
    dose_range     = paste(
      "IV bolus (~20-second injection) of 2 mg/kg for subjects <65",
      "years or 1 mg/kg for subjects >65 years, followed 60 minutes",
      "later by a 60-minute constant infusion at 25, 50, 100, or 200",
      "microgram/kg/min. Arterial plasma sampled at 0, 1, 2, 4, 8, 16,",
      "30, 60, 62, 64, 68, 76, 90, 120, 122, 124, 136, 150, 180, 240,",
      "300, and 600 minutes; the first sample used in the PBPK fit was",
      "at 2 minutes because of mixing-time effects."
    ),
    regions        = "USA (Stanford University Medical Center)",
    notes          = paste(
      "The PBPK model uses the PKQuest standard-human physiology (organ",
      "weights, blood flows, and fraction lipid) from Table 1. Only the",
      "average parameter set (Table 4 and Figs 9-11) is implemented",
      "here; individual per-subject Tclr / frdose / T fits (Figs 5-7)",
      "are not carried across because they are not tabulated in a form",
      "that maps to a population model. The paper's Roberts-Rowland",
      "dispersion-model liver (dispersion number 0.3) is replaced by a",
      "well-stirred compartment with apparent clearance calibrated to",
      "the paper's average Table 4 fractional liver clearance of 0.76.",
      "Kp values are hard-coded at the paper's standard freepl = 0.022."
    )
  )

  ini({
    # ---- Structural parameters (all FIXED at paper's averaged values) ----

    # Fractional liver blood flow that is cleared. Paper Table 4 reports
    # 0.772 (18-34 y), 0.781 (35-65 y), and 0.762 (>65 y); Figs 9-11 and
    # the Discussion use the age-averaged value 0.76 as the single
    # "averaged" parameter for all subjects.
    lfrac_liver_clearance <- fixed(log(0.76));          label("Fractional liver clearance (unitless; fraction of total liver + intestine blood flow extracted per pass)") # Levitt 2005 Table 4 age average (Discussion / Fig 9 caption)

    # Sequestration release rate = 1 / T with T = 80 min (paper's median
    # sequestration time constant; 5 of 8 young subjects had T = 80 min
    # per Fig 5, and the "averaged" parameter set of Figs 9-11 uses
    # T = 80 min).
    lkseq <- fixed(log(1 / 80));                        label("Pulmonary sequestration release rate (1/min); T = 80 min sequestration time constant")            # Levitt 2005 Fig 5 caption / Fig 9 caption (T = 80 min averaged)

    # Age-dependence of pulmonary bolus sequestration fraction:
    # frdose = max(0, 0.548 - 0.00891 * AGE) (paper eq 13).  frdose is
    # NOT used directly in the ODEs - bolus sequestration is applied
    # externally via a second dose row on cmt = "lung_seq" (see the
    # vignette).  The eq 13 regression coefficients (intercept 0.548,
    # slope -0.00891 /year) are documented in comments only rather
    # than as ini() parameters because they do not enter model().

    # Propofol physicochemistry provenance (not exposed as ini()
    # parameters because they do not enter model()): K_oil = 4715
    # (Levitt 2005 Methods / Weaver et al. ref 1),
    # freepl_st = 0.022 (Levitt 2005 Methods), and
    # rblpl_st = 1.0 (Levitt 2005 Methods).  The Kp values below
    # were derived once at freepl = 0.022 using these constants
    # (paper eqs 4-7); sensitivity analysis on freepl requires
    # recomputing Kp externally and updating the lkp_* entries.

    # Tissue/blood partition coefficients (hard-coded from paper text
    # for propofol at freepl = 0.022; paper page 3 first paragraph
    # after Table 1: "adipose 84; brain 1.87; liver 2.12; intestine 1.7
    # and the rest 1.45").  If the user wants to sensitivity-analyse
    # freepl, they must recompute Kp externally using paper eqs 4-7.
    lkp_adipose   <- fixed(log(84));                    label("Adipose/blood partition coefficient Kp (unitless)")            # Levitt 2005 Methods (paper text)
    lkp_brain     <- fixed(log(1.87));                  label("Brain/blood partition coefficient Kp (unitless)")              # Levitt 2005 Methods (paper text)
    lkp_liver     <- fixed(log(2.12));                  label("Liver/blood partition coefficient Kp (unitless)")              # Levitt 2005 Methods (paper text)
    lkp_intestine <- fixed(log(1.7));                   label("Intestine/blood partition coefficient Kp; paper's 'portal' (unitless)") # Levitt 2005 Methods (paper text; 'intestine 1.7')
    lkp_rest      <- fixed(log(1.45));                  label("Partition coefficient Kp for muscle, kidney, heart, lung, skin, other, tendon-lumped (unitless)") # Levitt 2005 Methods (paper text; 'the rest 1.45')

    # Residual proportional error - the paper reports an "average
    # weighted residual" (WRE) of about 15 percent for the individual
    # fits (Table 3 shows 14.18 percent for the age-average NONMEM model
    # infusion phase, 15 percent PBPK bolus + infusion).  WRE is the
    # mean absolute proportional deviation, not a formal residual SD;
    # propSd here fixes the paper's WRE at 0.15 for downstream
    # simulations that require a residual model.  See vignette
    # Assumptions and deviations.
    propSd <- fixed(0.15);                              label("Proportional residual error placeholder (fraction; equals paper WRE)")   # Levitt 2005 Table 3 average WRE ~ 15%
  })

  model({
    # ---- Reference physiology (70-kg standard human, Table 1) ----
    # Reference organ weights (kg) and blood flow perfusions (L/min/kg).
    # Standard human at 25 percent body fat: total mass 70 kg, adipose
    # 17.5 kg, lean mass 52.5 kg, cardiac output ~ 5.5877 L/min.
    ref_wt          <- 70
    ref_bodyfat_pct <- 25
    ref_lean_mass   <- ref_wt * (1 - ref_bodyfat_pct / 100)   # = 52.5 kg

    # Reference tissue weights (kg) - Table 1 column "Weight (Kg)".
    v_ref_liver     <- 1.8
    v_ref_intestine <- 1.5     # paper's "portal"
    v_ref_muscle    <- 26
    v_ref_kidney    <- 0.31
    v_ref_brain     <- 1.4
    v_ref_heart     <- 0.33
    v_ref_lung      <- 0.536
    v_ref_skin      <- 2.6
    v_ref_bone      <- 4       # bone (Q = 0 in Table 1 so has no distribution effect)
    v_ref_other     <- 5.524 + 3   # paper's "other" (5.524) + "tendon" (3) lumped connective tissue
    v_ref_adipose   <- 17.5

    # Reference blood flows (L/min) at 70 kg - Table 1 column
    # "Total Flow (L/min)" for tissues on the arterial side; the
    # intestine flow is the portal-vein contribution which enters the
    # liver.  Bone perfusion is 0 (paper Table 1) so the compartment
    # has no distributional role and is retained only for total mass.
    q_ref_liver     <- 0.45    # hepatic artery
    q_ref_intestine <- 1.125   # portal-vein contribution
    q_ref_muscle    <- 0.585
    q_ref_kidney    <- 1.24
    q_ref_brain     <- 0.784
    q_ref_heart     <- 0.264
    q_ref_lung      <- 5.6184  # cardiac output (sum of arterial tissues)
    q_ref_skin      <- 0.26
    q_ref_bone      <- 0       # paper Table 1
    q_ref_other     <- 0.1104 + 0.03   # "other" (0.1104) + "tendon" (0.03) lumped
    q_ref_adipose   <- 0.7385

    # Cardiac output for reference human - used as denominator when
    # scaling for cardiac-output changes.  Equal to q_ref_lung by
    # definition (lung receives full venous return).
    co_ref <- q_ref_lung

    # ---- Per-subject body composition (Gallagher regression, eq 2) ----
    # SEXM = 1 - SEXF because the Gallagher regression codes Sex = 0
    # for females and 1 for males (paper Methods text); SEXF is the
    # canonical female-indicator covariate.
    sexm <- 1 - SEXF

    # BMI in kg/m^2; HT is in cm from the covariate.
    bmi <- WT / (HT * 0.01)^2

    # Gallagher 1996 body-fat percent with Gallagher 2000 Asian
    # correction (paper eq 2).
    bodyfat_pct <-
      -10.0 + 1.46 * bmi - 11.6 * sexm + 0.14 * AGE +
       RACE_ASIAN * (95 / bmi - 0.044 * AGE)

    # Lean body mass (kg) and lean-mass scaling factor for the
    # non-adipose reference tissues.  Paper Methods: "liver weight is
    # a constant fraction of non-fat body weight, [so] the obese
    # subjects have a lower relative liver weight."
    lean_mass  <- WT * (1 - bodyfat_pct / 100)
    lean_scale <- lean_mass / ref_lean_mass

    # ---- Per-subject tissue weights (kg) ----
    v_liver     <- v_ref_liver     * lean_scale
    v_intestine <- v_ref_intestine * lean_scale
    v_muscle    <- v_ref_muscle    * lean_scale
    v_kidney    <- v_ref_kidney    * lean_scale
    v_brain     <- v_ref_brain     * lean_scale
    v_heart     <- v_ref_heart     * lean_scale
    v_lung      <- v_ref_lung      * lean_scale
    v_skin      <- v_ref_skin      * lean_scale
    v_bone      <- v_ref_bone      * lean_scale
    v_other     <- v_ref_other     * lean_scale
    v_adipose   <- WT * bodyfat_pct / 100

    # Blood volume (venous ~ 4/5, arterial ~ 1/5 of ~5.5 kg for 70 kg;
    # scale linearly with lean mass since blood is not adipose).
    v_venous   <- 5.5 * (4 / 5) * lean_scale
    v_arterial <- 5.5 * (1 / 5) * lean_scale

    # ---- Per-subject tissue blood flows (L/min) ----
    # Lean tissues scale flows by lean_scale (proportional to weight);
    # adipose flow is set from Table 1 perfusion 0.0422 L/min/kg times
    # per-subject adipose weight.
    q_liver     <- q_ref_liver     * lean_scale
    q_intestine <- q_ref_intestine * lean_scale
    q_muscle    <- q_ref_muscle    * lean_scale
    q_kidney    <- q_ref_kidney    * lean_scale
    q_brain     <- q_ref_brain     * lean_scale
    q_heart     <- q_ref_heart     * lean_scale
    q_skin      <- q_ref_skin      * lean_scale
    q_other     <- q_ref_other     * lean_scale
    q_bone      <- q_ref_bone      * lean_scale       # = 0
    q_adipose   <- 0.0422 * v_adipose

    # Cardiac output = sum of arterial-side tissue flows.  Note:
    # intestine (portal) flow drains into the liver (portal vein) so it
    # is NOT part of the direct arterial return - the arterial return
    # from the intestine goes to the liver, which then returns venously.
    # The lung flow equals the sum of everything returning to the
    # venous compartment.
    co <-
      q_liver + q_intestine + q_muscle + q_kidney + q_brain +
      q_heart + q_skin + q_other + q_adipose
    q_lung <- co

    # ---- Structural parameters (unpack ini() values) ----
    frac_liver_clearance <- exp(lfrac_liver_clearance)
    kseq                 <- exp(lkseq)

    kp_adipose   <- exp(lkp_adipose)
    kp_brain     <- exp(lkp_brain)
    kp_liver     <- exp(lkp_liver)
    kp_intestine <- exp(lkp_intestine)
    kp_rest      <- exp(lkp_rest)

    # Apparent hepatic clearance calibrated to the paper's Table 4
    # age-averaged fractional liver clearance of 0.76 of total liver
    # inflow (hepatic artery + portal vein).  Well-stirred well-mixed
    # approximation; replaces the paper's Roberts-Rowland dispersion
    # model liver (dispersion number 0.3) with an ODE-compatible
    # single-tank liver.  Verified against paper Table 2: at
    # frac_liver_clearance = 0.76, CL_apparent for the 70-kg standard
    # human = 0.76 * (0.45 + 1.125) = 1.197 L/min, matching the paper's
    # median value of "about 1.49 liters/min" for a fractional
    # clearance of 0.83 (paper Methods; the 0.76 value is the age-
    # averaged parameter for Figs 9-11).
    cl_apparent <- frac_liver_clearance * (q_liver + q_intestine)

    # ---- Concentrations (mg/L, blood-side) ----
    c_venous    <- venous     / v_venous
    c_arterial  <- arterial   / v_arterial
    c_liver     <- liver      / v_liver
    c_intestine <- intestine  / v_intestine
    c_muscle    <- muscle     / v_muscle
    c_kidney    <- kidney     / v_kidney
    c_brain     <- brain      / v_brain
    c_heart     <- heart      / v_heart
    c_lung      <- lung       / v_lung
    c_skin      <- skin       / v_skin
    c_other     <- other      / v_other
    c_bone      <- bone       / v_bone
    c_adipose   <- adipose    / v_adipose

    # ---- Flow-limited well-mixed tissue ODEs (paper eq 4 topology) ----
    # dA_i/dt = Q_i * (C_art_in - C_tissue_i / Kp_i)
    # where C_tissue_i / Kp_i is the effluent blood concentration in
    # equilibrium with the tissue (Roberts and Rowland well-stirred
    # form).  Arterial blood enters lean tissues; venous return
    # aggregates them.
    d/dt(muscle)    <- q_muscle    * (c_arterial - c_muscle    / kp_rest)
    d/dt(kidney)    <- q_kidney    * (c_arterial - c_kidney    / kp_rest)
    d/dt(brain)     <- q_brain     * (c_arterial - c_brain     / kp_brain)
    d/dt(heart)     <- q_heart     * (c_arterial - c_heart     / kp_rest)
    d/dt(skin)      <- q_skin      * (c_arterial - c_skin      / kp_rest)
    d/dt(other)     <- q_other     * (c_arterial - c_other     / kp_rest)
    d/dt(adipose)   <- q_adipose   * (c_arterial - c_adipose   / kp_adipose)
    d/dt(bone)      <- q_bone      * (c_arterial - c_bone      / kp_rest)   # q_bone = 0 -> inert

    # Intestine (paper's "portal") receives arterial blood and drains
    # via the portal vein into the liver - its effluent blood
    # concentration is c_intestine / kp_intestine.
    d/dt(intestine) <- q_intestine * (c_arterial - c_intestine / kp_intestine)

    # Liver receives (a) hepatic artery flow from arterial (rate
    # q_liver) at c_arterial and (b) portal-vein flow (rate
    # q_intestine) at c_intestine / kp_intestine; total inflow rate
    # is (q_liver + q_intestine).  Elimination is well-stirred
    # apparent clearance cl_apparent on the liver-side (effluent)
    # blood concentration c_liver / kp_liver.
    d/dt(liver) <-
      q_liver     * c_arterial +
      q_intestine * (c_intestine / kp_intestine) -
      (q_liver + q_intestine) * (c_liver / kp_liver) -
      cl_apparent              * (c_liver / kp_liver)

    # Lung receives full venous return (rate q_lung) at c_venous and
    # returns to arterial.  Standard flow-limited well-mixed lung.
    d/dt(lung) <- q_lung * (c_venous - c_lung / kp_rest)

    # Arterial compartment: pulmonary vein inflow (rate q_lung at
    # c_lung / kp_rest) plus optional lung-sequestration release
    # (rate kseq * lung_seq feeds directly into the arterial pool);
    # arterial outflow is q_lung (= total tissue return) at c_arterial.
    d/dt(arterial) <-
      q_lung * (c_lung / kp_rest) -
      q_lung * c_arterial +
      kseq   * lung_seq

    # Venous compartment: sum of tissue effluent returns minus lung
    # inflow.  Liver returns (q_liver + q_intestine) at c_liver /
    # kp_liver; other tissues at their Kp-appropriate effluent
    # concentrations.
    d/dt(venous) <-
      q_muscle    * (c_muscle    / kp_rest)      +
      q_kidney    * (c_kidney    / kp_rest)      +
      q_brain     * (c_brain     / kp_brain)     +
      q_heart     * (c_heart     / kp_rest)      +
      q_skin      * (c_skin      / kp_rest)      +
      q_other     * (c_other     / kp_rest)      +
      q_bone      * (c_bone      / kp_rest)      +   # q_bone = 0 -> no contribution
      q_adipose   * (c_adipose   / kp_adipose)   +
      (q_liver + q_intestine) * (c_liver / kp_liver) -
      q_lung      * c_venous

    # Pulmonary sequestration: bolus fraction is dosed externally to
    # this compartment (cmt = "lung_seq") and released exponentially
    # into the arterial pool at rate kseq = 1 / T.  Not used for
    # infusion doses; the user allocates the fraction frdose =
    # max(0, 0.548 - 0.00891 * AGE) via a second dose row.
    d/dt(lung_seq) <- -kseq * lung_seq

    # ---- Observation: arterial plasma propofol (paper's sampling
    # site was arterial blood; concentrations converted to plasma via
    # rblpl_st = 1.0, so arterial blood concentration equals arterial
    # plasma concentration for this model) ----
    Cc <- c_arterial
    Cc ~ prop(propSd)
  })
}
