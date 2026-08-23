Jiang_2025_filgrastim <- function() {
  description <- paste(
    "Semi-mechanistic population PK-PD model for filgrastim (recombinant",
    "human G-CSF) and peripheral-blood CD34+ cell mobilisation in healthy",
    "Korean male subjects given a single subcutaneous dose of 5 or",
    "10 ug/kg. PK is one-compartment with linear elimination preceded by a",
    "single transit compartment describing absorption delay: dose enters",
    "'depot' and moves at ktr_pk to the absorption compartment 'transit1',",
    "which is absorbed at ka into 'central'. Theory-based allometry scales",
    "CL (exponent 0.75) and V (exponent 1) by body weight. CD34+",
    "mobilisation is a modified Friberg model assuming continual entry of",
    "proliferating bone-marrow stem cells into peripheral blood through a",
    "single transit compartment: 'precursor1' (proliferating pool) ->",
    "'precursor2' (transit) -> 'circ' (circulating pool). A single rate",
    "constant ktr_pd serves simultaneously as the self-replication rate",
    "kprol, the transit rate and the circulating-pool elimination rate kel",
    "(constrained equal for identifiability), so the chain sits at",
    "steady state at the baseline count N0. G-CSF stimulates proliferation",
    "through an Emax function of plasma concentration, and the circulating",
    "count exerts negative feedback (N0/Circ)^gamma. Inter-occasion",
    "variability on ktr_pk and N0 is multiplexed by an OCC indicator over",
    "the two crossover periods."
  )
  reference <- paste(
    "Jiang X, Cha JS, Jin BH, Kim CO, Chae D. (2025). Population",
    "Pharmacokinetic-Pharmacodynamic Modeling of Granulocyte",
    "Colony-Stimulating Factor to Optimize Dosing and Timing for CD34+",
    "Cell Harvesting. Clin Transl Sci 18(1):e70121.",
    "doi:10.1111/cts.70121.",
    sep = " "
  )
  vignette <- "Jiang_2025_filgrastim"

  units <- list(
    time          = "h",
    dosing        = "ug",
    concentration = "ng/mL",
    CD34          = "cells/uL",
    notes         = paste(
      "Filgrastim dose in ug enters compartment 'depot' (the paper's dosing",
      "was 5 or 10 ug/kg subcutaneously, i.e. about 350 or 700 ug at the",
      "70 kg reference weight). Drug amounts are in ug and V is in L, so",
      "central / vc is in ug/L = ng/mL, matching the EC50 unit reported in",
      "Table 2. Bioavailability is not reported: only subcutaneous data were",
      "fitted, so CL and V are apparent values (CL/F and V/F) and F is",
      "implicitly 1. The three CD34+ chain states (precursor1, precursor2,",
      "circ) are cell counts in cells/uL, not amounts; CD34 = circ."
    )
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix.
  compartmentData <- list(
    depot       = list(analyte = "filgrastim", units = "ug", specimen = "administration site", verified = TRUE),
    transit1    = list(analyte = "filgrastim", units = "ug", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "filgrastim", units = "ug", specimen = "plasma", verified = TRUE),
    precursor1  = list(analyte = "CD34+ cells", units = "cells/uL", specimen = "not applicable", verified = TRUE),
    precursor2  = list(analyte = "CD34+ cells", units = "cells/uL", specimen = "not applicable", verified = TRUE),
    circ        = list(analyte = "CD34+ cells", units = "cells/uL", specimen = "not applicable", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "The only covariate retained in the final model (Jiang 2025",
        "Results 3.4). Theory-based allometry on CL (exponent 0.75 FIX) and",
        "V (exponent 1 FIX), Formulas (4) and (5) and Table 2. IMPORTANT:",
        "Formulas (4) and (5) are printed WITHOUT a normalising reference",
        "weight (CL_i = theta_CL,pop x WT^0.75; V_i = theta_V,pop x WT).",
        "Taken literally those give CL = 29.7 L/h and V = 293 L at 70 kg,",
        "which would put Cmax near 2 ng/mL after a 10 ug/kg dose, whereas",
        "the observed Cmax in Figure 3a is 35-75 ng/mL. The Table 2 values",
        "are therefore typical values at a reference weight and the",
        "relationships are encoded in the standard weight-normalised form",
        "(WT / 70)^exponent; 70 kg is the theory-based-allometry standard",
        "(reference 28, Anderson and Holford) and also matches the study",
        "median weight of 70.40 kg (Table 1). See the vignette Errata.",
        "Observed range 60.80-85.50 kg (Table 1)."
      ),
      source_name        = "WT"
    ),
    OCC = list(
      description        = "Integer-valued crossover-period indicator for inter-occasion variability multiplexing.",
      units              = "(count)",
      type               = "count",
      reference_category = NULL,
      notes              = paste(
        "Values 1 and 2 identify the treatment period within subject. The",
        "study was a randomised, open-label, two-way crossover single-dose",
        "bioequivalence trial in which each subject received Neupogen in one",
        "period and the biosimilar Leucostim in the other (Jiang 2025",
        "Methods 2.1), and inter-occasion variability was included",
        "'due to the crossover study design' (Methods 2.3). Decomposed",
        "inside model() into binary indicators oc1 and oc2 that multiplex the",
        "IOV etas on log-ktr_pk and log-N0 (Table 2 IOV rows). The paper",
        "found no significant PK or PD difference between the two",
        "formulations and pooled them, so OCC carries period only and is not",
        "a formulation effect."
      ),
      source_name        = "OCC"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = "Screened in the stepwise covariate analysis (Jiang 2025 Methods 2.3, Table 1) but not retained; body weight was the only significant covariate. Observed range 20-44 years."
    ),
    HT = list(
      description = "Body height",
      units       = "cm",
      type        = "continuous",
      notes       = "Screened but not retained (Jiang 2025 Methods 2.3, Table 1). Observed range 163.60-189.00 cm."
    ),
    BMI = list(
      description = "Body mass index",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = "Screened but not retained (Jiang 2025 Methods 2.3, Table 1). Observed range 19.00-24.90 kg/m^2 (an inclusion criterion restricted BMI to 18.5-25.0)."
    ),
    WBC = list(
      description = "White blood cell count",
      units       = "10^3 cells/uL",
      type        = "continuous",
      notes       = "Screened but not retained (Jiang 2025 Methods 2.3, Table 1). Observed range 3.62-9.98 x 10^3/uL."
    ),
    RBC = list(
      description = "Red blood cell count",
      units       = "10^6 cells/uL",
      type        = "continuous",
      notes       = "Screened but not retained (Jiang 2025 Methods 2.3, Table 1). Observed range 4.53-5.94 x 10^6/uL."
    ),
    HGB = list(
      description = "Hemoglobin",
      units       = "g/dL",
      type        = "continuous",
      notes       = "Screened but not retained (Jiang 2025 Methods 2.3, Table 1). Observed range 12.60-17.30 g/dL."
    ),
    HCT = list(
      description = "Hematocrit",
      units       = "%",
      type        = "continuous",
      notes       = "Screened but not retained (Jiang 2025 Methods 2.3, Table 1). Observed range 39.50-50.10%."
    ),
    NEUT = list(
      description = "Absolute neutrophil count at baseline",
      units       = "10^3 cells/uL",
      type        = "continuous",
      notes       = "Screened but not retained (Jiang 2025 Methods 2.3, Table 1, reported as ANC). Observed range 2.20-5.87 x 10^3/uL; an inclusion criterion required 2-7 x 10^3/uL."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 53L,
    n_studies      = 1L,
    age_range      = "20-44 years (median 31; Table 1)",
    age_median     = "31 years",
    weight_range   = "60.80-85.50 kg (median 70.40; Table 1)",
    weight_median  = "70.40 kg",
    sex_female_pct = 0,
    race_ethnicity = "Korean",
    disease_state  = paste(
      "Healthy adult male volunteers; no history of recombinant G-CSF",
      "administration, baseline absolute neutrophil count within",
      "2-7 x 10^3/uL, body weight at least 60 kg and BMI 18.5-25.0 kg/m^2."
    ),
    dose_range     = paste(
      "Single subcutaneous dose of 5 ug/kg (Part A, n = 26) or 10 ug/kg",
      "(Part B, n = 27). Each subject received Neupogen in one period and",
      "the biosimilar Leucostim in the other; the two formulations were",
      "pooled."
    ),
    regions        = "Republic of Korea (Severance Hospital, Seoul)",
    notes          = paste(
      "Trial registration NCT02725086; a randomised, open-label, two-way",
      "crossover, single-dose Phase I bioequivalence study of Neupogen",
      "(Amgen) versus the biosimilar Leucostim (Dong-A Pharmaceutical).",
      "1378 plasma filgrastim concentrations and 982 CD34+ cell counts.",
      "PK sampling at 0, 1, 2, 3, 4, 5, 6, 8, 12, 16, 24, 36 and 48 h;",
      "CD34+ sampling at 0, 24, 48, 72, 96, 120, 168, 240 and 312 h.",
      "Filgrastim measured by ELISA (Quantikine Human G-CSF, R&D Systems),",
      "calibration 0-625 pg/mL with LLOQ 5.05 pg/mL; CD34+ counts by flow",
      "cytometry (BD FACSVerse with the BD Stem Cell Enumeration Kit).",
      "Estimation software Monolix Suite 2024R1 (SAEM); 1000-replicate",
      "bootstrap and 500-replicate VPC. In addition to the covariates",
      "recorded in covariatesDataExcluded the stepwise screen also tested",
      "platelet count and the leukocyte differential percentages",
      "(neutrophil, lymphocyte, monocyte, eosinophil and basophil percent);",
      "none was retained, and none has a canonical covariate-register name,",
      "so they are recorded here rather than as covariatesDataExcluded",
      "entries. Because only two dose levels were studied and linear",
      "elimination sufficed, the authors explicitly declined to extrapolate",
      "to higher doses (Discussion)."
    )
  )

  ini({
    # -----------------------------------------------------------------
    # PK structural parameters. CL and V are apparent (CL/F, V/F): only
    # subcutaneous data were fitted and no bioavailability term is
    # reported. Both are typical values at the 70 kg reference weight --
    # see the covariateData WT note and the vignette Errata for why the
    # weight relationship is encoded in normalised form.
    # -----------------------------------------------------------------
    lvc       <- log(4.19)                ; label("Log V (apparent central volume of distribution) at WT 70 kg (L)")                            # Jiang 2025 Table 2 (V = 4.19 L, RSE 3.64%; bootstrap median 3.87, 95% CI 2.417-5.383)
    lcl       <- log(1.23)                ; label("Log CL (apparent clearance) at WT 70 kg (L/h)")                                              # Jiang 2025 Table 2 (CL = 1.23 L/h, RSE 3.34%; bootstrap median 1.24, 95% CI 1.156-1.337)
    lka       <- log(0.56)                ; label("Log ka (first-order absorption rate from the absorption compartment) (1/h)")                  # Jiang 2025 Table 2 (ka = 0.56 1/h, RSE 5.00%; bootstrap median 0.56, 95% CI 0.389-0.767)
    lktr      <- log(0.25)                ; label("Log ktr_pk (transit rate constant of G-CSF from depot to the absorption compartment) (1/h)")  # Jiang 2025 Table 2 (ktr_pk = 0.25 1/h, RSE 4.52%; bootstrap median 0.23, 95% CI 0.180-0.328)

    e_wt_cl   <- fixed(0.75)              ; label("Allometric exponent on (WT / 70 kg) for CL (unitless)")                                      # Jiang 2025 Table 2 (betaCL_logWT = 0.75 FIX) and Formula (4); theory-based allometry per reference 28
    e_wt_vc   <- fixed(1)                 ; label("Allometric exponent on (WT / 70 kg) for V (unitless)")                                       # Jiang 2025 Table 2 (betaV_logWT = 1 FIX) and Formula (5); theory-based allometry per reference 28

    # -----------------------------------------------------------------
    # CD34+ mobilisation (modified Friberg). Figure 1 annotates
    # kprol(=ktr_pd) and kel(=ktr_pd): the self-replication rate of the
    # proliferating pool and the elimination rate of the circulating pool
    # are both constrained equal to the transit rate ktr_pd, so only one
    # PD rate constant is estimated. Discussion: "The transit and
    # elimination rate constants were assumed to be identical (estimated
    # at 0.059/h) due to parameter identifiability issues." With all three
    # rates equal the chain is at steady state when every state equals
    # N0, which is what sets the initial conditions in model().
    # -----------------------------------------------------------------
    lrbase    <- log(6.19)                ; label("Log N0 (baseline circulating CD34+ cell count) (cells/uL)")                                  # Jiang 2025 Table 2 (N0 = 6.19, RSE 6.82%; bootstrap median 5.72, 95% CI 3.491-7.804) -- but see the vignette Errata: the paper's own Figures 2e/2h and 3c/3d imply a typical baseline near 1.3-2 cells/uL
    lktr_cd34 <- log(0.059)               ; label("Log ktr_pd (CD34+ transit rate; also kprol and kel) (1/h)")                                  # Jiang 2025 Table 2 (ktr_pd = 0.059 1/h, RSE 3.79%; bootstrap median 0.073, 95% CI 0.055-0.124) and Figure 1 annotations kprol(=ktr_pd), kel(=ktr_pd)
    lemax     <- log(1.15)                ; label("Log Emax (maximum fractional stimulation of CD34+ proliferation) (unitless)")                 # Jiang 2025 Table 2 (Emax = 1.15, RSE 8.11%; bootstrap median 0.6, 95% CI 0.274-1.308)
    lec50     <- log(0.11)                ; label("Log EC50 (plasma G-CSF concentration giving half-maximal stimulation) (ng/mL)")               # Jiang 2025 Table 2 (EC50 = 0.11 ng/mL, RSE 32.2%; bootstrap median 0.0081, 95% CI 0.000164-0.197 -- the least precisely estimated parameter)
    lgamma    <- log(0.2)                 ; label("Log gamma (exponent of the (N0 / Circ) negative-feedback term) (unitless)")                   # Jiang 2025 Table 2 (Gamma = 0.2, RSE 8.42%; bootstrap median 0.16, 95% CI 0.0884-0.234) and Figure 1 (Feedback = (Circ0 / Circ(t))^gamma)

    # -----------------------------------------------------------------
    # Inter-individual variability. Table 2 reports the IIV rows as
    # "SD (CV)" (footnote b), i.e. the first number is the log-scale SD
    # omega and the parenthesised number is the corresponding
    # CV = sqrt(exp(omega^2) - 1) in percent. ini() takes variances, so
    # each entry below is log(CV^2 + 1); the CV column is used because it
    # carries more significant figures than the two-significant-digit SD
    # column (e.g. omega_CL 0.23 <-> CV 23.31% reproduces exactly, while
    # omega_gamma 0.34 rounds up from 0.3388 <-> CV 34.89%).
    # No IIV was estimated on ka or EC50.
    # -----------------------------------------------------------------
    etalvc       ~ log(0.1089^2 + 1)  # Jiang 2025 Table 2: omega V      = 0.11 (CV 10.89%), RSE 35.5%
    etalcl       ~ log(0.2331^2 + 1)  # Jiang 2025 Table 2: omega CL     = 0.23 (CV 23.31%), RSE 10.9%
    etalktr      ~ log(0.1288^2 + 1)  # Jiang 2025 Table 2: omega ktr_pk = 0.13 (CV 12.88%), RSE 39.5%
    etalrbase    ~ log(0.4181^2 + 1)  # Jiang 2025 Table 2: omega N0     = 0.4  (CV 41.81%), RSE 11.5%
    etalktr_cd34 ~ log(0.1646^2 + 1)  # Jiang 2025 Table 2: omega ktr_pd = 0.16 (CV 16.46%), RSE 26.7%
    etalemax     ~ log(0.2461^2 + 1)  # Jiang 2025 Table 2: omega Emax   = 0.24 (CV 24.61%), RSE 16.0%
    etalgamma    ~ log(0.3489^2 + 1)  # Jiang 2025 Table 2: omega gamma  = 0.34 (CV 34.89%), RSE 21.7%

    # -----------------------------------------------------------------
    # Inter-occasion variability across the two crossover periods, on
    # ktr_pk and N0 only (Table 2 IOV rows, same SD (CV) convention). One
    # eta per occasion sharing a single variance, multiplexed by OCC in
    # model(); occasion 2 is fixed to the occasion-1 variance because the
    # paper reports a single IOV magnitude per parameter.
    # -----------------------------------------------------------------
    etaiov_ktr_1   ~ log(0.2194^2 + 1)         # Jiang 2025 Table 2: gamma ktr_pk = 0.22 (CV 21.94%), RSE 14.3%
    etaiov_ktr_2   ~ fixed(log(0.2194^2 + 1))  # occasion 2 shares the occasion-1 IOV variance
    etaiov_rbase_1 ~ log(0.1495^2 + 1)         # Jiang 2025 Table 2: gamma N0      = 0.15 (CV 14.95%), RSE 14.6%
    etaiov_rbase_2 ~ fixed(log(0.1495^2 + 1))  # occasion 2 shares the occasion-1 IOV variance

    # -----------------------------------------------------------------
    # Residual error. PK: combined additive and proportional (Results
    # 3.3, Formula 3). PD: proportional (Results 3.3, Formula 2).
    # -----------------------------------------------------------------
    addSd        <- 0.4                   ; label("Additive residual SD on plasma G-CSF concentration Cc (ng/mL)")                              # Jiang 2025 Table 2 (magnitude of additive error for G-CSF = 0.4, RSE 4.27%)
    propSd       <- 0.32                  ; label("Proportional residual SD on plasma G-CSF concentration Cc")                                  # Jiang 2025 Table 2 (magnitude of proportional error for G-CSF = 0.32, RSE 2.87%)
    propSd_CD34  <- 0.35                  ; label("Proportional residual SD on the circulating CD34+ cell count")                                # Jiang 2025 Table 2 (magnitude of proportional error for CD34+ = 0.35, RSE 2.79%)
  })

  model({
    # -----------------------------------------------------------------
    # 1. Occasion indicators and IOV multiplexing (two crossover periods)
    # -----------------------------------------------------------------
    oc1       <- (OCC == 1)
    oc2       <- (OCC == 2)
    iov_ktr   <- oc1 * etaiov_ktr_1   + oc2 * etaiov_ktr_2
    iov_rbase <- oc1 * etaiov_rbase_1 + oc2 * etaiov_rbase_2

    # -----------------------------------------------------------------
    # 2. Individual parameters. Formulas (4) and (5) in weight-normalised
    #    form at the 70 kg theory-based-allometry reference.
    # -----------------------------------------------------------------
    vc        <- exp(lvc + etalvc) * (WT / 70) ^ e_wt_vc
    cl        <- exp(lcl + etalcl) * (WT / 70) ^ e_wt_cl
    ka        <- exp(lka)
    ktr       <- exp(lktr + etalktr + iov_ktr)

    rbase     <- exp(lrbase + etalrbase + iov_rbase)
    ktr_cd34  <- exp(lktr_cd34 + etalktr_cd34)
    emax      <- exp(lemax + etalemax)
    ec50      <- exp(lec50)
    gamma     <- exp(lgamma + etalgamma)

    # -----------------------------------------------------------------
    # 3. Micro-constant
    # -----------------------------------------------------------------
    kel       <- cl / vc

    # -----------------------------------------------------------------
    # 4. Plasma concentration, drug effect and negative feedback.
    #    Defined before the ODE block because the CD34+ equations consume
    #    them. Figure 1: E = Emax * Cp / (EC50 + Cp) stimulates kprol, and
    #    Feedback = (Circ0 / Circ(t))^gamma damps it.
    # -----------------------------------------------------------------
    Cc        <- central / vc
    edrug     <- emax * Cc / (ec50 + Cc)
    feedback  <- (rbase / circ) ^ gamma

    # -----------------------------------------------------------------
    # 5. ODE system.
    #    depot      : subcutaneous filgrastim, the paper's dosing site (ug)
    #    transit1   : the paper's "Absorption Compartment" (ug)
    #    central    : plasma filgrastim (ug)
    #    precursor1 : proliferating bone-marrow stem-cell pool (cells/uL)
    #    precursor2 : the paper's "Transit 2" mobilisation delay (cells/uL)
    #    circ       : circulating peripheral-blood CD34+ pool (cells/uL)
    #    kprol and the circulating-pool kel are both ktr_cd34 (Figure 1).
    # -----------------------------------------------------------------
    d/dt(depot)      <- -ktr * depot
    d/dt(transit1)   <-  ktr * depot - ka * transit1
    d/dt(central)    <-  ka  * transit1 - kel * central

    d/dt(precursor1) <-  ktr_cd34 * precursor1 * (1 + edrug) * feedback - ktr_cd34 * precursor1
    d/dt(precursor2) <-  ktr_cd34 * precursor1 - ktr_cd34 * precursor2
    d/dt(circ)       <-  ktr_cd34 * precursor2 - ktr_cd34 * circ

    # -----------------------------------------------------------------
    # 6. Initial conditions. With kprol = ktr_pd = kel the CD34+ chain is
    #    at steady state when all three states equal the baseline count.
    # -----------------------------------------------------------------
    precursor1(0)    <- rbase
    precursor2(0)    <- rbase
    circ(0)          <- rbase

    # -----------------------------------------------------------------
    # 7. Observations
    # -----------------------------------------------------------------
    Cc   ~ add(addSd) + prop(propSd)
    CD34 <- circ
    CD34 ~ prop(propSd_CD34)
  })
}
