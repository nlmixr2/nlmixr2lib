Zeng_2026_atezolizumab <- function() {
  description <- "Two-compartment population PK model for intravenous atezolizumab (anti-PD-L1 IgG1) in patients with metastatic non-small cell lung cancer, with sigmoidal time-varying clearance (Zeng 2026): every fixed effect, variance and residual-error term is held FIXED for simulation, transcribed from the FDA CDER clinical pharmacology and biopharmaceutics review for BLA 761041Orig1s000. The source publication is a three-way software benchmark (NONMEM, RxODE, Pumas) whose supplement carries the complete NONMEM control stream plus the equivalent RxODE and Pumas model objects."
  reference <- paste(
    "Zeng Y, Arisa O, Corvalan N, Bateman F, Schmidt K, Peer C, Figg WD.",
    "A comparison of pharmacometric software programs for atezolizumab",
    "population pharmacokinetic simulation.",
    "Eur J Clin Pharmacol. 2026;82:42. doi:10.1007/s00228-025-03974-9.",
    "Parameter values duplicated by the authors from:",
    "US Food and Drug Administration, Center for Drug Evaluation and Research.",
    "Clinical Pharmacology and Biopharmaceutics Review,",
    "application number 761041Orig1s000 (atezolizumab, TECENTRIQ), 2016.",
    "The same model was used for the dose re-optimization analysis in",
    "Peer CJ, Schmidt KT, Arisa O, et al. J Clin Pharmacol. 2023;63(6):672-680.",
    sep = " "
  )
  vignette <- "Zeng_2026_atezolizumab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. verified = TRUE: the source reports simulated *serum*
  # atezolizumab concentrations (Zeng 2026 Figs. 1-3 captions, Table 2 in
  # ug/mL) and the control stream scales the central compartment by V1
  # (S1 = V1), so both states hold atezolizumab as an amount in mg.
  compartmentData <- list(
    central     = list(analyte = "atezolizumab", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "atezolizumab", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power scaling on CL and Vc with reference weight 77 kg (Zeng 2026 Supplement, NONMEM $PK: TVCL includes (WT/77)**THETA(7) and TVV1 includes (WT/77)**THETA(11)). The virtual cohort draws WT from an age- and sex-dependent normal distribution (Zeng 2026 Supplement, population-generation code).",
      source_name        = "WT"
    ),
    ALB = list(
      description        = "Baseline serum albumin",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power scaling on CL and Vc with reference 40 g/L (Zeng 2026 Supplement, NONMEM $PK: (ALB/40)**THETA(8) on CL and (ALB/40)**THETA(12) on V1). Reported in g/L (SI): Zeng 2026 Table 1 gives the virtual-cohort median as 42.0-42.2 g/L and the supplement draws it as rnorm(n, 42, sd = 3.5). The Methods sentence 'a mean of 42 g/dL' is a unit typo in the source text; Table 1 and the covariate reference of 40 both establish g/L.",
      source_name        = "ALB"
    ),
    TUMSZ = list(
      description        = "Baseline tumor size (sum of diameters of target lesions)",
      units              = "mm",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power scaling on CL with reference 63 mm (Zeng 2026 Supplement, NONMEM $PK: (BTS/63)**THETA(9)). Zeng 2026 Table 1 reports the virtual-cohort baseline tumor size in mm (median 67.5-69.3 mm) and the supplement draws it as rlnorm(n, 4.2, sd = 0.7), i.e. a median of exp(4.2) = 66.7 mm. The Methods sentence 'baseline tumor burden ... 4.2 mm^3' conflates the log-scale meanlog with a volume; Table 1 and the generating code establish mm.",
      source_name        = "BTS"
    ),
    SEXF = list(
      description        = "Biological sex indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Multiplicative factor on Vc (0.896) and Vp (0.707) for females relative to males (Zeng 2026 Supplement, NONMEM $PK: (THETA(13)**ISF) on V1 and (THETA(14)**ISF) on V2). The control-stream column is named ISF and the population-generation code sets SEX = 1 for female, 0 for male, so ISF is an is-female indicator and maps directly onto the canonical SEXF with no value transformation.",
      source_name        = "ISF"
    ),
    ADA_POS = list(
      description        = "Antidrug-antibody status",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (ADA-negative)",
      notes              = "Multiplicative factor of 1.158 on CL for ADA-positive patients (Zeng 2026 Supplement, NONMEM $PK: (THETA(10)**ATAG)). Renamed from the source column ATAG (anti-therapy antibody) to the canonical ADA_POS. Drawn as rbinom(n, 1, 0.4) in the virtual cohort, giving the 38.3-43.3 percent ADA-positive fractions in Zeng 2026 Table 1.",
      source_name        = "ATAG"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Baseline age",
      units       = "years",
      type        = "continuous",
      notes       = "Not a covariate of the PK model. Age enters only the virtual-population generator, where it drives body weight (Zeng 2026 Supplement: BW = 65 + 0.75 * (AGE - 40) + N(0, 3.5) for females and 85 + 0.75 * (AGE - 40) + N(0, 10) for males, with AGE ~ U(20, 80)). Documented here so the cohort-construction provenance is preserved without declaring an unused model covariate."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 1000,
    n_studies      = 1,
    age_range      = "20-80 years",
    age_median     = "44.6-51 years across the three software implementations",
    weight_range   = "43.2-136 kg",
    weight_median  = "79.6-83.2 kg across the three software implementations",
    sex_female_pct = 51.2,
    disease_state  = "Metastatic non-small cell lung cancer (the underlying population PK model was developed in NSCLC patients per the FDA CDER review of BLA 761041Orig1s000).",
    dose_range     = "840, 1200 or 1680 mg IV; standard 1200 mg q3w and three extended-interval regimens (2 loading cycles at 840 mg q2w, 1200 mg q3w or 1680 mg q4w followed by 840 mg q6w), 7 cycles total",
    notes          = "The 1000 subjects are a SIMULATED virtual population, not an estimation dataset: this model was not fitted to observed data in the source publication. Covariates were generated in R, Julia and NONMEM from the distributions given in the Zeng 2026 Supplement (seed 12345): AGE ~ U(20, 80) rounded; SEX Bernoulli(0.5) with 1 = female; WT normal about an age- and sex-dependent mean; ADA_POS ~ Bernoulli(0.4); ALB ~ N(42, 3.5) g/L; TUMSZ ~ lognormal(meanlog 4.2, sdlog 0.7) mm. Baseline albumin and tumor-size distributions were anchored on the medians reported by Stroh et al. Clin Pharmacol Ther 2017;102(2):305-312. Realised cohort demographics are in Zeng 2026 Table 1. The size of the FDA estimation dataset is not reported in this publication."
  )

  ini({
    # EVERY parameter below is FIXED. The source publication performs
    # simulation only: the NONMEM control stream in the Zeng 2026 Supplement
    # carries FIX on all 17 $THETA entries, all 4 $OMEGA entries and $SIGMA,
    # and the values are transcribed from the FDA CDER clinical pharmacology
    # and biopharmaceutics review for BLA 761041Orig1s000. None was estimated
    # from data in this paper.
    #
    # Time is in DAYS. The NONMEM and RxODE implementations in the supplement
    # are on a day scale; the Pumas implementation restates the same model on
    # an hour scale (tvcl 0.009583 L/h * 24 = 0.23 L/day; tvq 0.025125 L/h *
    # 24 = 0.603 L/day; tvT50 1507.1999 h / 24 = 62.8 day), which confirms the
    # day-scale reading.

    # Structural parameters -- typical values at WT 77 kg, ALB 40 g/L,
    # TUMSZ 63 mm, male, ADA-negative.
    lcl <- fixed(log(0.23));  label("Clearance at reference covariates (CL, L/day)")             # Zeng 2026 Supplement, NONMEM $THETA(3)
    lvc <- fixed(log(3.25));  label("Central volume of distribution at reference (Vc, L)")       # Zeng 2026 Supplement, NONMEM $THETA(4)
    lq  <- fixed(log(0.603)); label("Intercompartmental clearance (Q, L/day)")                   # Zeng 2026 Supplement, NONMEM $THETA(5)
    lvp <- fixed(log(2.88));  label("Peripheral volume of distribution at reference (Vp, L)")    # Zeng 2026 Supplement, NONMEM $THETA(6)

    # Sigmoidal time-varying clearance: cl <- cl_base * exp(cl_time_max_i *
    # t^gamma / (t50^gamma + t^gamma)). cl_time_max is negative, so clearance
    # falls with time on treatment, by a factor of exp(-0.193) = 0.824 as
    # t >> cl_t50. Kept on the natural scale because it is signed.
    cl_time_max   <- fixed(-0.193); label("Maximal fractional change in log-CL over time (unitless)")     # Zeng 2026 Supplement, NONMEM $THETA(15) IMAX
    cl_t50   <- fixed(62.8);   label("Time at which half of cl_time_max is reached (day)")           # Zeng 2026 Supplement, NONMEM $THETA(16) T50
    cl_time_hill <- fixed(2.67);   label("Hill / sigmoidicity exponent of time on CL (unitless)")        # Zeng 2026 Supplement, NONMEM $THETA(17) GAMMA

    # Covariate effects on CL. WT, ALB and TUMSZ enter as power terms; ADA
    # status enters as a multiplicative factor raised to the 0/1 indicator,
    # exactly as the control stream parameterizes it.
    e_wt_cl    <- fixed(0.668);  label("Power exponent of WT on CL (unitless)")                          # Zeng 2026 Supplement, NONMEM $THETA(7)
    e_alb_cl   <- fixed(-0.901); label("Power exponent of ALB on CL (unitless)")                         # Zeng 2026 Supplement, NONMEM $THETA(8)
    e_tumsz_cl <- fixed(0.116);  label("Power exponent of TUMSZ on CL (unitless)")                       # Zeng 2026 Supplement, NONMEM $THETA(9)
    e_ada_cl   <- fixed(1.158);  label("Multiplicative factor of ADA-positive status on CL (unitless)")  # Zeng 2026 Supplement, NONMEM $THETA(10)

    # Covariate effects on the volumes.
    e_wt_vc   <- fixed(0.533);  label("Power exponent of WT on Vc (unitless)")                           # Zeng 2026 Supplement, NONMEM $THETA(11)
    e_alb_vc  <- fixed(-0.345); label("Power exponent of ALB on Vc (unitless)")                          # Zeng 2026 Supplement, NONMEM $THETA(12)
    e_sexf_vc <- fixed(0.896);  label("Multiplicative factor of female sex on Vc (unitless)")            # Zeng 2026 Supplement, NONMEM $THETA(13)
    e_sexf_vp <- fixed(0.707);  label("Multiplicative factor of female sex on Vp (unitless)")            # Zeng 2026 Supplement, NONMEM $THETA(14)

    # IIV, all FIXED. The four $OMEGA entries are variances on the log scale
    # and are exact squares of the standard deviations the FDA review reports:
    # 0.263^2 = 0.069169 (CL), 0.172^2 = 0.029584 (Vc), 0.352^2 = 0.123904
    # (Vp), 0.897^2 = 0.804609 (cl_time_max).
    etalcl         ~ fixed(0.069169)  # Zeng 2026 Supplement, NONMEM $OMEGA ETA(1), IIV on CL
    etalvc         ~ fixed(0.029584)  # Zeng 2026 Supplement, NONMEM $OMEGA ETA(2), IIV on V1
    etalvp         ~ fixed(0.123904)  # Zeng 2026 Supplement, NONMEM $OMEGA ETA(3), IIV on V2
    etacl_time_max ~ fixed(0.804609)  # Zeng 2026 Supplement, NONMEM $OMEGA ETA(4), IIV on IMAX

    # Combined additive plus proportional residual error, FIXED. The control
    # stream computes W = SQRT(THETA(1)**2 * IPRED**2 + THETA(2)**2) with
    # $SIGMA 1 FIX, i.e. THETA(1) and THETA(2) are the proportional and
    # additive standard deviations.
    propSd <- fixed(0.034); label("Proportional residual error (fraction)")   # Zeng 2026 Supplement, NONMEM $THETA(1) with $ERROR W
    addSd  <- fixed(18.1);  label("Additive residual error (ug/mL)")          # Zeng 2026 Supplement, NONMEM $THETA(2) with $ERROR W
  })
  model({
    # Individual parameters (Zeng 2026 Supplement, NONMEM $PK). Reference
    # covariates: WT 77 kg, ALB 40 g/L, TUMSZ 63 mm, male, ADA-negative.
    cl_base <- exp(lcl + etalcl) *
      (WT    / 77)^e_wt_cl *
      (ALB   / 40)^e_alb_cl *
      (TUMSZ / 63)^e_tumsz_cl *
      e_ada_cl^ADA_POS

    vc <- exp(lvc + etalvc) *
      (WT  / 77)^e_wt_vc *
      (ALB / 40)^e_alb_vc *
      e_sexf_vc^SEXF

    vp <- exp(lvp + etalvp) * e_sexf_vp^SEXF
    q  <- exp(lq)

    # Time-varying clearance. IIV enters multiplicatively on the signed
    # cl_time_max, as IMAX = THETA(15) * EXP(ETA(4)) in the control stream,
    # so the sign of the effect is preserved for every subject.
    cl_time_max_i <- cl_time_max * exp(etacl_time_max)
    cl <- cl_base * exp(cl_time_max_i * t^cl_time_hill /
                          (cl_t50^cl_time_hill + t^cl_time_hill))

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                  k12 * central - k21 * peripheral1

    # Observation: dose in mg, volume in L -> mg/L = ug/mL (S1 = V1 in the
    # control stream).
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
