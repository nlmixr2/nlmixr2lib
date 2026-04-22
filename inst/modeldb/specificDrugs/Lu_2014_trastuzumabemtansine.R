Lu_2014_trastuzumabemtansine <- function() {
  description <- "Linear two-compartment population PK model of trastuzumab emtansine (T-DM1, anti-HER2 antibody-drug conjugate) with first-order elimination from the central compartment in patients with HER2-positive locally advanced or metastatic breast cancer (Lu 2014)"
  reference <- "Lu D, Girish S, Gao Y, Wang B, Yi J-H, Guardino E, Samant M, Cobleigh M, Rimawi M, Conte P, Jin JY. Population pharmacokinetics of trastuzumab emtansine (T-DM1), a HER2-targeted antibody-drug conjugate, in patients with HER2-positive metastatic breast cancer: clinical implications of the effect of covariates. Cancer Chemother Pharmacol. 2014;74(2):399-410. doi:10.1007/s00280-014-2500-2"
  vignette <- "Lu_2014_trastuzumabemtansine"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline value. Power effect on CL (exponent 0.49) and Vc (exponent 0.596); reference 70 kg per Lu 2014 Eq. 5.",
      source_name        = "weight"
    ),
    HER2_ECD = list(
      description        = "Baseline serum HER2 shed extracellular domain concentration",
      units              = "ng/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline value. Power effect on CL (exponent 0.035); reference 25 ng/mL per Lu 2014 Eq. 5. Source column 'ECD' maps to the canonical scope-specific HER2_ECD covariate.",
      source_name        = "ECD"
    ),
    ALB = list(
      description        = "Baseline serum albumin concentration",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline value (SI units; not g/dL). Power effect on CL (exponent -0.423); reference 41 g/L per Lu 2014 Eq. 5. Source column 'ALBU' maps to the canonical ALB covariate.",
      source_name        = "ALBU"
    ),
    TUMSZ = list(
      description        = "Baseline sum of the longest dimension of target lesions (RECIST)",
      units              = "mm",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline value. Power effect on CL (exponent 0.052); reference 90 mm per Lu 2014 Eq. 5 (source reference 9 cm * 10 = 90 mm). Source column 'TMBD' is in cm; convert values to mm on data ingestion (TUMSZ_mm = TMBD_cm * 10) so the canonical unit and per-model reference agree.",
      source_name        = "TMBD"
    ),
    TRAST_BL = list(
      description        = "Baseline serum trastuzumab concentration from prior trastuzumab-containing therapy",
      units              = "ug/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline value. Enters the CL equation linearly (not log-transformed) via exp(coef * TRAST_BL); reference value 0 ug/mL corresponds to trastuzumab-naive patients. Coefficient -0.002 ug/mL^-1 per Lu 2014 Eq. 5. Source column 'TBL' maps to the canonical scope-specific TRAST_BL covariate. Lu 2014 observed TBL = 0 ug/mL at the 5th percentile and 54 ug/mL at the 95th percentile.",
      source_name        = "TBL"
    ),
    AST = list(
      description        = "Baseline serum aspartate aminotransferase activity",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline value. Power effect on CL (exponent 0.071); reference 27 U/L per Lu 2014 Eq. 5. Source column 'AST' maps to the canonical AST covariate.",
      source_name        = "AST"
    )
  )

  population <- list(
    n_subjects     = 671L,
    n_studies      = 5L,
    phase_mix      = "1 phase I (TDM3569g), 3 phase II (TDM4258g, TDM4374g, TDM4450g), and 1 phase III (TDM4370g / EMILIA) studies.",
    regimen_mix    = "643 (95.8%) patients on 3.6 mg/kg IV every 3 weeks; 28 (4.2%) patients on weekly (qw) regimens.",
    n_observations = 9934L,
    disease_state  = "HER2-positive locally advanced or metastatic breast cancer (MBC).",
    dose_range     = "Clinical doses 2.4-4.8 mg/kg IV q3w (labeled regimen 3.6 mg/kg q3w); phase I included lower doses (>=0.3 mg/kg) where nonlinear clearance was observed at <=1.2 mg/kg.",
    weight_range   = "5th-95th percentile 49-98 kg (per Lu 2014 Table 2); mean Asian 60.5 kg vs non-Asian 71.6 kg.",
    age_range      = "Adults with MBC; sensitivity subgroups <65, 65-75, >75 years (per Lu 2014 Supplemental Table 4).",
    sex_female_pct = NA_real_,
    race_ethnicity = "Global enrollment including Asian and non-Asian patients (Lu 2014 Supplemental Table 3). Race and region were tested but were not statistically significant covariates.",
    regions        = "Global / multi-regional.",
    external_validation = "Phase II study TDM4688g (n = 51; 505 concentration-time data points) used for external validation.",
    reference_subject  = "70 kg, HER2_ECD 25 ng/mL, ALB 41 g/L, TUMSZ 90 mm (9 cm), TRAST_BL 0 ug/mL, AST 27 U/L (per Lu 2014 Table 2 footnote a and Eq. 5).",
    notes          = "Baseline demographic and clinical characteristics are reported in Lu 2014 Supplemental Table 3. Observations below the minimum quantifiable concentration (0.04-0.06 ug/mL) were excluded (7.27% of data points). Age, race, geographic region, and renal function were tested but were not statistically significant covariates on CL or Vc."
  )

  ini({
    # Structural parameters (typical values for the reference subject defined
    # in Lu 2014 Table 2 footnote a). Paper parameterization uses clearances on
    # an hour basis via exp(theta)*24 -> L/day, reported in Table 1; we work
    # directly in L/day so the stored typical value is log(CL_Lperday).
    lcl  <- log(0.676);  label("Clearance at reference covariates (CL, L/day)")                      # Lu 2014 Table 1 (exp(theta1)*24 = 0.676 L/day)
    lvc  <- log(3.127);  label("Central volume of distribution at reference (Vc, L)")               # Lu 2014 Table 1 (exp(theta2) = 3.127 L)
    lq   <- log(1.534);  label("Intercompartmental clearance (Q, L/day)")                           # Lu 2014 Table 1 (exp(theta3)*24 = 1.534 L/day)
    lvp  <- log(0.66);   label("Peripheral volume of distribution (Vp, L)")                         # Lu 2014 Table 1 (exp(theta4) = 0.66 L)

    # Covariate-effect parameters (Lu 2014 Eq. 5 and Table 1).
    # Continuous covariates enter as power models normalized to the Table 2
    # reference values (70 kg, ECD 25 ng/mL, ALB 41 g/L, TMBD 9 cm = 90 mm,
    # AST 27 U/L). TBL enters linearly on the log-CL scale with reference 0
    # (no prior trastuzumab).
    e_wt_cl    <-  0.49;   label("Power exponent of WT on CL (unitless)")                           # Lu 2014 Table 1, theta6
    e_wt_vc    <-  0.596;  label("Power exponent of WT on Vc (unitless)")                           # Lu 2014 Table 1, theta5
    e_ecd_cl   <-  0.035;  label("Power exponent of HER2_ECD on CL (unitless)")                     # Lu 2014 Table 1, theta7
    e_alb_cl   <- -0.423;  label("Power exponent of ALB on CL (unitless)")                          # Lu 2014 Table 1, theta8
    e_tumsz_cl <-  0.052;  label("Power exponent of TUMSZ on CL (unitless)")                        # Lu 2014 Table 1, theta9
    e_trast_cl <- -0.002;  label("Linear coefficient of TRAST_BL on log-CL (per ug/mL)")            # Lu 2014 Table 1, theta10
    e_ast_cl   <-  0.071;  label("Power exponent of AST on CL (unitless)")                          # Lu 2014 Table 1, theta11

    # Inter-individual variability. Lu 2014 Table 1 reports the NONMEM IIV as
    # %CV on log-normal parameters with Var(eta) = log(CV^2 + 1). The CL-Vc
    # pair is a correlated block (covariance 0.011 on the log scale; implied
    # correlation ~0.500). IIV on Q and Vp is retained as reported despite the
    # large shrinkage (49.7% and 36.1%) noted in the Results section.
    etalcl + etalvc ~ c(0.03585,
                        0.011, 0.01351)  # CL 19.11% CV, Vc 11.66% CV, cov 0.011 -- Lu 2014 Table 1
    etalq  ~ 1.4513                       # Q 180.8% CV -- Lu 2014 Table 1
    etalvp ~ 0.4415                       # Vp 74.50% CV -- Lu 2014 Table 1

    # Residual error. Lu 2014 Table 1 reports a single proportional residual
    # error of 31.56% CV (Sigma row). propSd is stored as the SD on the
    # proportional scale.
    propSd <- 0.3156; label("Proportional residual error (SD, fraction)")                           # Lu 2014 Table 1 (Sigma = 31.56% CV)
  })
  model({
    # Individual PK parameters with covariate adjustments (Lu 2014 Eq. 5).
    # Reference covariate values: WT 70 kg, HER2_ECD 25 ng/mL, ALB 41 g/L,
    # TUMSZ 90 mm (9 cm), TRAST_BL 0 ug/mL, AST 27 U/L. Power effects act on
    # log-transformed continuous covariates in NONMEM, which translate to
    # (COV / ref)^exponent in linear space. TBL is linear on the log-CL scale,
    # i.e. exp(e_trast_cl * TRAST_BL) in linear space.
    cl <- exp(lcl + etalcl) *
      (WT       /   70)^e_wt_cl *
      (HER2_ECD /   25)^e_ecd_cl *
      (ALB      /   41)^e_alb_cl *
      (TUMSZ    /   90)^e_tumsz_cl *
      (AST      /   27)^e_ast_cl *
      exp(e_trast_cl * TRAST_BL)

    vc <- exp(lvc + etalvc) *
      (WT / 70)^e_wt_vc

    q  <- exp(lq + etalq)
    vp <- exp(lvp + etalvp)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1

    # Concentration: dose in mg, volume in L -> mg/L = ug/mL.
    Cc <- central / vc

    Cc ~ prop(propSd)
  })
}
