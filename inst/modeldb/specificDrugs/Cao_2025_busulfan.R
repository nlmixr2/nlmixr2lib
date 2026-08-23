Cao_2025_busulfan <- function() {
  description <- "Semi-mechanistic two-compartment IV popPK model for busulfan in pediatric allogeneic hematopoietic cell transplantation recipients, with normal-fat-mass allometric scaling on CL/Q and Vc/Vp, a postmenstrual-age sigmoid maturation function on CL, and a baseline-normalised glutathione pool that makes elimination time-varying, scaled by baseline glutathione S-transferase enzyme activity (Cao 2025)."
  reference <- "Cao D, Qian X, Wang P, Zheng X, Huang S, Wei Z, Jiang W, Yu L, Jiang X, Yu Y, Mao J, Zhai X. Semi-mechanistic population pharmacokinetic model incorporating glutathione S-transferase activity for personalized busulfan dosing in pediatric allogeneic hematopoietic cell transplantation. Front Pharmacol. 2025;16:1632588. doi:10.3389/fphar.2025.1632588"
  vignette <- "Cao_2025_busulfan"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against the Supplementary Text S2 NONMEM control
  # stream ($MODEL NCOMP=3: CENTRAL / PERIPH / GSH; S1 = V1/1000 confirms the
  # central amount is in mg while the fitted concentration scale was ng/mL).
  compartmentData <- list(
    central     = list(analyte = "busulfan", units = "mg", specimen = "whole blood", verified = TRUE),
    peripheral1 = list(analyte = "busulfan", units = "mg", specimen = "whole blood", verified = TRUE),
    gsh_pool    = list(
      analyte  = "glutathione",
      units    = "(relative to baseline; dimensionless, initialised at 1)",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed at baseline. Combined with FFM to build parameter-specific normal fat mass, NFM = FFM + Ffat * (WT - FFM) (Supplementary Text S1 Eq. 4). The adult reference is a 70 kg subject, whose FFM is 56.1 kg (see FFM notes).",
      source_name        = "WT"
    ),
    FFM = list(
      description        = "Fat-free mass (Janmahasatian 2005 semi-mechanistic model, derived from body weight, height, and sex)",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed at baseline. Supplementary Text S1 Eqs. 2-3 give FFM(male) = 9270 * WT / (6680 + 216 * BMI) and FFM(female) = 9270 * WT / (8780 + 244 * BMI), with BMI in kg/m^2. The control stream hardcodes the adult reference FFM as 56.1 kg, which is exactly the male formula evaluated at WT = 70 kg and HT = 1.76 m; that identity is what fixes the reference and the sex-specific formula used for the standard subject.",
      source_name        = "FFM"
    ),
    AGE = list(
      description        = "Postnatal age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Combined with GA to form postmenstrual age inside model(): PMA (weeks) = AGE * 365 / 7 + GA, exactly as written in the Supplementary Text S2 $PK block. Note the source uses 365/7 = 52.14 weeks per year, not 365.25/7.",
      source_name        = "AGE"
    ),
    GA = list(
      description        = "Gestational age at birth",
      units              = "weeks",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Second component of postmenstrual age. Table 3 reports per-subject gestational ages of 37.1-39.1 weeks for the ten virtual-trial subjects; the cohort was term-born.",
      source_name        = "GAGE"
    ),
    GST_BL_NMOL_MIN_ML = list(
      description        = "Baseline whole-blood glutathione S-transferase enzyme activity, measured once per subject by micro-quartz colorimetry (GST-catalysed conjugation of glutathione to 1-chloro-2,4-dinitrobenzene, read at 340 nm)",
      units              = "nmol/min/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed; sampled once, before the first busulfan dose (Methods 2.2). Enters as a power function on the glutathione-depletion scaling factor, sdep_gsh * (GST_BL_NMOL_MIN_ML / 9.2)^e_gst_sdep_gsh, where 9.2 nmol/min/mL is the model-development cohort median (Table 1). The reference value equalling the reported cohort median is a transcription check on the centering constant. Observed range 0.9-20.7 nmol/min/mL in the development set and 3.7-12.4 in the evaluation set.",
      source_name        = "GST"
    ),
    OCC = list(
      description        = "Dosing-occasion index for inter-occasion variability on clearance",
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = "Four IOV occasions. The source control stream keys them off the dose number (OCC1 = dose 1, OCC2 = dose 5, OCC3 = dose 11, OCC4 = dose 12), which are precisely the doses around which busulfan was sampled (Methods 2.1). This model file re-indexes them to the canonical consecutive OCC = 1, 2, 3, 4 so the column follows the register; records outside a sampled occasion should carry OCC = 0, which zeroes every indicator and leaves CL at its between-subject value.",
      source_name        = "OCC"
    )
  )

  # Screened in the stepwise covariate search (Supplementary Table S3) but not
  # retained in the final model, so they are documentation only and are
  # deliberately absent from model().
  covariatesDataExcluded <- list(
    HCT = list(
      description = "Hematocrit",
      units       = "%",
      type        = "continuous",
      notes       = "Screened on CL as linear / power / exponential; best dOFV -3.3, below the forward-inclusion threshold of 3.84 (Supplementary Table S3 rows 2-4)."
    ),
    TBIL = list(
      description = "Total bilirubin",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Screened on CL; best dOFV -0.3 (Supplementary Table S3 rows 5-7)."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened on CL; best dOFV -1.7 (Supplementary Table S3 rows 8-10)."
    ),
    ALB = list(
      description = "Serum albumin",
      units       = "g/L",
      type        = "continuous",
      notes       = "Screened on CL; best dOFV -3.3, below the 3.84 threshold (Supplementary Table S3 rows 11-13)."
    ),
    CRP = list(
      description = "C-reactive protein",
      units       = "mg/L",
      type        = "continuous",
      notes       = "Screened on CL; dOFV 0.002, i.e. no effect (Supplementary Table S3 rows 14-16)."
    ),
    CRCL = list(
      description = "Estimated glomerular filtration rate",
      units       = "mL/min/1.73 m^2",
      type        = "continuous",
      notes       = "Screened on CL; best dOFV -0.75 (Supplementary Table S3 rows 17-19). Consistent with busulfan being almost entirely metabolised, with only about 2% excreted unchanged in urine.",
      source_name = "eGFR"
    ),
    CONMED_FLUDARABINE = list(
      description = "Concomitant fludarabine administration",
      units       = "(binary)",
      type        = "binary",
      notes       = "The only covariate to pass forward inclusion (dOFV -10.8, Supplementary Table S3 row 21), but it failed backward elimination: removing it raised the OFV by only 5.4, below the 10.83 retention threshold, so it was excluded from the final model (Results 3.2)."
    ),
    DIS_PID = list(
      description = "Primary immunodeficiency disease indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened on CL as a categorical effect; dOFV -1.1 (Supplementary Table S3 row 20).",
      source_name = "PID"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 65,
    n_studies      = 1,
    n_centers      = 1,
    n_observations = 636,
    age_range      = "0.2-14.1 years (model-development set; median 1.4 years, 21 subjects under 1 year)",
    weight_range   = "2.9-29.5 kg (model-development set; median 9.9 kg)",
    height_range   = "52.0-147.0 cm (model-development set; median 76.0 cm)",
    sex_female_pct = 27.7,
    disease_state  = "Pediatric recipients of allogeneic hematopoietic cell transplantation receiving intravenous busulfan as part of a myeloablative preparative regimen",
    dose_range     = "0.8-1.2 mg/kg per dose by weight band, given as a 2 h IV infusion every 6 h for 12 doses",
    regions        = "China (Children's Hospital of Fudan University, Shanghai)",
    gst_range      = "0.9-20.7 nmol/min/mL (model-development set; median 9.2)",
    notes          = "Prospective single-center cohort enrolled August 2020 to November 2021 (Methods 2.1, Table 1). 636 whole-blood busulfan concentrations from 65 subjects; 536 samples from 55 subjects were used for model development and 100 samples from 10 subjects were held out for external evaluation. One further participant was excluded for disease progression (Figure 1). Sampling was at 2, 2.5, 3, 4 and 6 h after dose 1, pre-dose before doses 6 and 12, and at 2/4/8 h (Group A, n = 33) or 2/6/12 h (Group B, n = 32) after dose 12."
  )

  ini({
    # Structural parameters, standardised to a 70 kg adult (Table 2; identical to
    # the Supplementary Text S2 $THETA block, records 1-4).
    lcl <- log(9.57);  label("Typical clearance at the adult reference size (CL, L/h)")                          # Table 2 / Text S2 $THETA(1)
    lvc <- log(28.2);  label("Typical central volume of distribution at the adult reference size (Vc, L)")       # Table 2 / Text S2 $THETA(2)
    lq  <- log(8.16);  label("Typical intercompartmental clearance at the adult reference size (Q, L/h)")        # Table 2 / Text S2 $THETA(3)
    lvp <- log(16.1);  label("Typical peripheral volume of distribution at the adult reference size (Vp, L)")    # Table 2 / Text S2 $THETA(4)

    # Fat-mass fractions defining parameter-specific normal fat mass. Estimated
    # here (unlike Lawson 2022, which fixed them a priori from McCune 2014).
    ffat_cl <- 0.905;  label("Fraction of fat mass contributing to NFM on CL and Q (unitless)")                  # Table 2 / Text S2 $THETA(5)
    ffat_vc <- 0.687;  label("Fraction of fat mass contributing to NFM on Vc and Vp (unitless)")                 # Table 2 / Text S2 $THETA(6)

    # Sigmoid maturation of CL on postmenstrual age (Equation 1).
    tm50_mat <- 45.0;  label("Postmenstrual age at which CL maturation reaches 50% of the adult value (weeks)")  # Table 2 / Text S2 $THETA(7)
    hill_mat <- 1.11;  label("Hill coefficient of the CL maturation function (unitless)")                        # Table 2 / Text S2 $THETA(8)

    # Glutathione-pool coupling. sdep_gsh scales busulfan metabolic flux to the
    # change in the normalised glutathione pool; it was not estimated here but
    # carried over from Langenhorst 2020.
    sdep_gsh <- fixed(0.00259)
    label("Scaling factor coupling busulfan metabolism to the glutathione pool (L/mg)")                          # Table 2 / Text S2 $THETA(9), value from Langenhorst 2020
    e_gst_sdep_gsh <- 0.28
    label("Power exponent on baseline GST activity (relative to 9.2 nmol/min/mL) for sdep_gsh (unitless)")       # Table 2 / Text S2 $THETA(10)

    # Between-subject variability. Table 2 reports these as percentages that are
    # sqrt(omega) * 100, so omega = (pct/100)^2. There is deliberately no BSV on
    # Q (Methods 2.3: "BSV modeled via log-normal distributions for all
    # parameters, except Q").
    etalcl ~ 0.0539   # Table 2 BSV CL 23.2% -> 0.232^2 = 0.0539
    etalvc ~ 0.0244   # Table 2 BSV Vc 15.6% -> 0.156^2 = 0.0244; matches Text S2 $OMEGA record 2 exactly
    etalvp ~ 0.16     # Table 2 BSV Vp 40.0% -> 0.400^2 = 0.16;   matches Text S2 $OMEGA record 3 exactly

    # Inter-occasion variability on CL over four sampled occasions. The source
    # encodes this as $OMEGA BLOCK(1) 0.0115 followed by three SAME blocks, i.e.
    # one shared magnitude across all four occasions.
    etaiov_cl_1 ~ 0.0115        # Table 2 IOV on CL 10.7% -> 0.107^2 = 0.0115; Text S2 $OMEGA BLOCK(1)
    etaiov_cl_2 ~ fixed(0.0115) # Text S2 $OMEGA BLOCK(1) SAME
    etaiov_cl_3 ~ fixed(0.0115) # Text S2 $OMEGA BLOCK(1) SAME
    etaiov_cl_4 ~ fixed(0.0115) # Text S2 $OMEGA BLOCK(1) SAME

    # Residual error, combined proportional + additive. Text S2 $SIGMA gives
    # 0.0124 and 276; sqrt(0.0124) = 0.1114 and sqrt(276) = 16.61, matching
    # Table 2's "Proportional 11.1%" and "Additional 16.6". The additive term is
    # on the fitted ng/mL scale (the control stream sets S1 = V1/1000), so it is
    # 16.61 ng/mL = 0.01661 mg/L in this model's declared mg/L units.
    #
    # Deviation: Text S2 $ERROR writes Y = F * EXP(EPS(1)) + EPS(2), i.e. a
    # log-normal rather than linear proportional term. rxode2 5.1.7 parses
    # `lnorm() + add()` but cannot solve it ("cannot find additive standard
    # deviation for 'Cc'"), and every model in this library must solve, so the
    # linear proportional form is used. This matches how the paper itself
    # describes the residual (Results 3.2, "combining proportional and additive
    # models provided the best results for the RUV"; Table 2 row label
    # "Proportional (%)"), and at this magnitude exp(eps) and 1 + eps differ by
    # under 1% of the one-sigma multiplier.
    propSd <- 0.1114;  label("Proportional residual standard deviation (fraction)")                              # Text S2 $SIGMA(1) = 0.0124
    addSd  <- 0.01661; label("Additive residual standard deviation (mg/L)")                                      # Text S2 $SIGMA(2) = 276 (ng/mL)^2
  })

  model({
    # ---- Body composition: parameter-specific normal fat mass -----------------
    # NFM = FFM + Ffat * (WT - FFM) (Supplementary Text S1 Eq. 4). Ffat = 0
    # collapses NFM to FFM; Ffat = 1 collapses it to total body weight. The
    # standard subject is a 70 kg adult whose FFM is 56.1 kg, so the reference
    # NFM is itself Ffat-dependent (Supplementary Text S2 $PK NFMSTD_CL /
    # NFMSTD_V).
    fat        <- WT - FFM
    nfm_cl     <- FFM + ffat_cl * fat
    nfm_v      <- FFM + ffat_vc * fat
    nfm_std_cl <- 56.1 + ffat_cl * (70 - 56.1)
    nfm_std_v  <- 56.1 + ffat_vc * (70 - 56.1)

    # Theory-based allometric exponents: 0.75 on the flow terms, 1 on the volumes.
    size_cl <- (nfm_cl / nfm_std_cl)^0.75
    size_v  <- nfm_v / nfm_std_v

    # ---- Maturation of clearance on postmenstrual age (Equation 1) ------------
    pma  <- AGE * 365 / 7 + GA
    fmat <- 1 / (1 + (pma / tm50_mat)^(-hill_mat))

    # ---- Inter-occasion variability on clearance ------------------------------
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    iov_cl <- oc1 * etaiov_cl_1 + oc2 * etaiov_cl_2 + oc3 * etaiov_cl_3 + oc4 * etaiov_cl_4

    # ---- Individual parameters ------------------------------------------------
    # Q carries the CL size factor (not the volume one) and no maturation term,
    # per Supplementary Text S2 $PK: Q = THETA(3) * FSIZCL.
    cl <- exp(lcl + etalcl + iov_cl) * size_cl * fmat
    vc <- exp(lvc + etalvc) * size_v
    q  <- exp(lq)           * size_cl
    vp <- exp(lvp + etalvp) * size_v

    # ---- Glutathione coupling (Equation 2) ------------------------------------
    # Baseline GST activity acts as a power function on the depletion scaling
    # factor, centered on the development-cohort median of 9.2 nmol/min/mL.
    sdep <- sdep_gsh * (GST_BL_NMOL_MIN_ML / 9.2)^e_gst_sdep_gsh

    kel  <- cl / vc
    k12  <- q  / vc
    k21  <- q  / vp
    kgsh <- sdep / vc

    # The glutathione pool is normalised to 1 at baseline. Busulfan elimination
    # is proportional to the pool, so the pool multiplies k10 in the central
    # equation and the busulfan metabolic flux drives the pool's own equation.
    # Both Equation 2 and the Supplementary Text S2 $DES block write this term
    # with a POSITIVE sign, so the pool grows and elimination accelerates over a
    # treatment course; see the vignette's Errata section, which reconciles this
    # against the paper's "GSH depletion" narrative and shows that the positive
    # sign is the one that reproduces the paper's own Table 3 simulations.
    gsh_pool(0) <- 1
    d/dt(central)     <- -(k12 + kel * gsh_pool) * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    d/dt(gsh_pool)    <-  kgsh * gsh_pool * kel * central

    Cc <- central / vc
    Cc ~ prop(propSd) + add(addSd)
  })
}
