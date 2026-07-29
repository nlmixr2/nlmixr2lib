Agerso_2004_desmopressin <- function() {
  description <- paste(
    "Three-compartment population PK model for intravenous desmopressin",
    "with simultaneous plasma and urinary-amount outputs. Systemic",
    "clearance is split into renal and non-renal components, each",
    "modulated linearly by creatinine clearance (CRCL), in healthy",
    "subjects and patients with varying degrees of renal impairment",
    "(Agerso 2004)."
  )
  reference <- paste(
    "Agerso H, Seiding Larsen L, Riis A, Lovgren U, Karlsson MO,",
    "Senderovitz T. Pharmacokinetics and renal excretion of desmopressin",
    "after intravenous administration to healthy subjects and renally",
    "impaired patients. Br J Clin Pharmacol. 2004;58(4):352-358.",
    "doi:10.1111/j.1365-2125.2004.02175.x."
  )
  vignette <- "Agerso_2004_desmopressin"
  paper_specific_etas <- c("etalvc_vp2")
  units <- list(time = "hour", dosing = "ng", concentration = "pg/mL")

  covariateData <- list(
    CRCL = list(
      description        = paste(
        "Creatinine clearance, BSA-normalised to 1.73 m^2.",
        "Agerso 2004 derived CRCL from a 24-h urine collection using",
        "(Cr_urine x UV / dt) / Cr_serum x 1.73 / BSA, where BSA is from",
        "the Gehan and George formula. Population CRCL range 10-144",
        "mL/min/1.73 m^2; reference value 60 mL/min/1.73 m^2 (population",
        "median)."
      ),
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed at study entry (single 24-h collection at baseline).",
        "Linear effect on both CL_r and CL_nr expressed as fractional",
        "change per mL/min from the population median of 60 mL/min/1.73",
        "m^2. The linear form is the paper's fitted relationship; it",
        "breaks down at CRCL < ~2.5 mL/min where the renal-clearance",
        "factor would become non-positive. Simulation outside the",
        "observed range (10-144 mL/min/1.73 m^2) should be capped."
      ),
      source_name        = "CrCL"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 22L,
    n_studies      = 1L,
    age_range      = "49-68 years",
    weight_range   = "BMI 23-32 kg/m^2; height 152-184 cm",
    sex_female_pct = NA_real_,
    disease_state  = paste(
      "Single-centre cohort stratified by renal function:",
      "group 1 normal (n=6, CrCL ~103 mL/min),",
      "group 2 mild impairment (n=5, ~72 mL/min),",
      "group 3 moderate impairment (n=5, ~37 mL/min),",
      "group 4 severe impairment (n=6, ~16 mL/min).",
      "None of the subjects were on dialysis."
    ),
    dose_range     = "Single 2 ug IV bolus (4 ug/mL, 0.5 mL over 1 min)",
    regions        = "Germany (single centre)",
    notes          = paste(
      "Pharmacokinetic analysis excluded 2 of 24 enrolled subjects whose",
      "plasma concentration-time profiles showed an absorption phase",
      "incompatible with IV administration (extravascular dosing",
      "suspected). Parameter table summarises N = 22. Plasma assay LLOQ",
      "2.5 pg/mL; urine assay LLOQ 20 pg/mL. See Methods 'Study design'",
      "and Methods 'Protocol' of Agerso 2004 for the design and the",
      "Results section / Table 2 for per-subject demographics."
    )
  )

  ini({
    # --- Structural disposition parameters (Table 1). Both clearance
    # components are reported at the population-median CRCL of 60
    # mL/min/1.73 m^2; the CRCL effect terms below apply the linear
    # fractional adjustment per mL/min deviation from 60.
    lcl_renal  <- log(2.41);  label("Renal clearance at CRCL=60 (L/h)")           # Table 1: CL_r = 2.41 L/h (RSE 7%)
    lcl_nonren <- log(3.91);  label("Non-renal clearance at CRCL=60 (L/h)")       # Table 1: CL_nr = 3.91 L/h (RSE 7%)
    lvc        <- log(13.9);  label("Central volume V1 (L)")                       # Table 1: V1 = 13.9 L (RSE 7%)
    lvp        <- log(12.5);  label("Peripheral volume V2 (L)")                    # Table 1: V2 = 12.5 L (RSE 8%)
    lvp2       <- log(11.7);  label("Second peripheral volume V3 (L)")             # Table 1: V3 = 11.7 L (RSE 8%)
    lq         <- log(4.43);  label("Intercompartmental clearance Q2 (L/h)")       # Table 1: Q2 = 4.43 L/h (RSE 21%)
    lq2        <- log(27.5);  label("Intercompartmental clearance Q3 (L/h)")       # Table 1: Q3 = 27.5 L/h (RSE 19%)

    # --- Covariate effects: CRCL on each clearance component. The paper
    # reports the relationship as a percentage change per 1 mL/min change
    # in CRCL evaluated at the population median CRCL = 60 mL/min.
    # Encoded here as a linear multiplicative model:
    #   CL_r  = CL_r_ref  * (1 + e_crcl_cl_renal  * (CRCL - 60))
    #   CL_nr = CL_nr_ref * (1 + e_crcl_cl_nonren * (CRCL - 60))
    e_crcl_cl_renal  <- 0.0174;   label("CRCL effect on CL_r (fraction per mL/min)")   # Table 1: 1.74% per 1 mL/min (RSE 2%)
    e_crcl_cl_nonren <- 0.00933;  label("CRCL effect on CL_nr (fraction per mL/min)")  # Table 1: 0.933% per 1 mL/min (RSE 12%)

    # --- IIV (log-normal). omega^2 = log(CV^2 + 1).
    # CL_r has a separate, uncorrelated eta (CV 24%, RSE 58%).
    # CL_nr and a shared V1/V3 eta are jointly distributed with
    # correlation 0.86 (Table 1 footnote ).
    #   var(etalcl_nonren) = log(1 + 0.30^2) = 0.08618    # Table 1 CV 30% (RSE 43%)
    #   var(etalvc_vp2)    = log(1 + 0.31^2) = 0.09175    # Table 1 CV 31% (RSE 45%, V1 and V3 rows)
    #   cov                = 0.86 * sqrt(0.08618 * 0.09175) = 0.07648
    etalcl_renal                 ~ 0.05599                            # CV 24%
    etalcl_nonren + etalvc_vp2 ~ c(0.08618,
                                     0.07648, 0.09175)                # CV 30%, r = 0.86, CV 31%

    # --- Residual error (combined additive + proportional), separate for
    # plasma concentrations and cumulative urinary amounts.
    propSd            <- 0.29;    label("Proportional residual SD on plasma Cc (fraction)")    # Table 1: prop error on Cp 29% (RSE 10%)
    addSd             <- 0.27;    label("Additive residual SD on plasma Cc (pg/mL)")            # Table 1: add error on Cp 0.27 pg/mL (RSE 59%)
    propSd_urineAmt   <- 0.47;    label("Proportional residual SD on urinary amount (fraction)") # Table 1: prop error on urine 47% (RSE 21%)
    addSd_urineAmt    <- 4140;    label("Additive residual SD on urinary amount (pg)")           # Table 1: add error on urine 4140 pg (RSE 21%)
  })

  model({
    # --- Individual disposition parameters. CRCL effects act on CL
    # components via a linear (1 + e * (CRCL - 60)) multiplier; the
    # remaining structural parameters are unaffected (the paper found no
    # other significant covariates).
    crcl_effect_renal  <- 1 + e_crcl_cl_renal  * (CRCL - 60)
    crcl_effect_nonren <- 1 + e_crcl_cl_nonren * (CRCL - 60)

    cl_renal  <- exp(lcl_renal  + etalcl_renal)  * crcl_effect_renal
    cl_nonren <- exp(lcl_nonren + etalcl_nonren) * crcl_effect_nonren
    vc        <- exp(lvc        + etalvc_vp2)                          # shared eta with vp2
    vp        <- exp(lvp)
    vp2       <- exp(lvp2       + etalvc_vp2)                          # shared eta with vc
    q         <- exp(lq)
    q2        <- exp(lq2)

    # --- Micro-constants for the three-compartment model.
    kel_r  <- cl_renal  / vc          # renal elimination from central
    kel_nr <- cl_nonren / vc          # non-renal elimination from central
    k12    <- q  / vc
    k21    <- q  / vp
    k13    <- q2 / vc
    k31    <- q2 / vp2

    # --- ODE system. The urine compartment accumulates the renal-cleared
    # amount; total elimination from central is the sum of renal and
    # non-renal components.
    d/dt(central)     <- -(kel_r + kel_nr) * central -
                          k12 * central + k21 * peripheral1 -
                          k13 * central + k31 * peripheral2
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    d/dt(peripheral2) <-  k13 * central - k31 * peripheral2
    d/dt(urine)       <-  kel_r * central

    # --- Observations.
    # Cc        : plasma concentration. Dose in ng / vc in L -> ng/L = pg/mL.
    # urineAmt  : cumulative amount in urine, reported in the paper in pg.
    #             The urine compartment accumulates in ng (dose units), so
    #             multiply by 1000 to match Table 1's additive error scale.
    Cc       <- central / vc
    urineAmt <- urine * 1000

    Cc       ~ add(addSd)          + prop(propSd)
    urineAmt ~ add(addSd_urineAmt) + prop(propSd_urineAmt)
  })
}
