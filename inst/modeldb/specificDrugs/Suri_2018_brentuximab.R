Suri_2018_brentuximab <- function() {
  description <- "Coupled population PK model for brentuximab vedotin antibody-drug conjugate (ADC) and its released payload monomethyl auristatin E (MMAE) in 380 patients with CD30-positive malignancies (Hodgkin lymphoma, systemic anaplastic large-cell lymphoma, mycosis fungoides, primary cutaneous ALCL) pooled from six clinical studies including the phase III ALCANZA study (Suri 2018). ADC is described by a linear 3-compartment model with zero-order input and first-order elimination; MMAE by a 2-compartment model with first-order elimination, fed from ADC by (a) a saturable target-binding flux Kd*Target*ADC (initial Target = 1, irreversibly depleted) and (b) a proteolytic flux FM*exp(-ALFM*tad)*K10*ADC whose conversion fraction declines as a function of time after the most recent dose. Both fluxes accumulate in an intermediate Lag compartment that empties to MMAE central with rate Klag (FM is fixed to 1)."
  reference <- "Suri A, Mould DR, Liu Y, Jang G, Venkatakrishnan K. Population PK and Exposure-Response Relationships for the Antibody-Drug Conjugate Brentuximab Vedotin in CTCL Patients in the Phase III ALCANZA Study. Clin Pharmacol Ther. 2018;104(5):989-999. doi:10.1002/cpt.1037. PMID 29377077."
  vignette <- "Suri_2018_brentuximab"
  units <- list(
    time          = "hour",
    dosing        = "umol",
    concentration = "umol/L"
  )

  covariateData <- list(
    BSA = list(
      description        = "Baseline body surface area",
      units              = "m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline. Power effects on ADC central volume V1 (exponent 1.27), ADC clearance CL (exponent 0.457), MMAE central volume V1 (exponent 0.89), and MMAE clearance CL (exponent 2.81). Reference 1.865 m^2 (Suri 2018 Table 1 overall-population mean; supplement 1 states continuous covariates were 'normalized for the population mean'). Pooled-population BSA range 1.264-2.858 m^2.",
      source_name        = "BSA"
    ),
    ALB = list(
      description        = "Baseline serum albumin concentration",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline. Power effects on ADC clearance CL (exponent -0.496) and MMAE clearance CL (exponent 0.982). Reference 36.81 g/L (Suri 2018 Table 1 overall-population mean). Pooled-population ALB range 17-53 g/L.",
      source_name        = "ALB"
    ),
    TBILI = list(
      description        = "Baseline total serum bilirubin concentration",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline. Power effect on MMAE clearance CL only (exponent -0.1). Reference 7.78 umol/L (Suri 2018 Table 1 overall-population mean). Pooled-population TBILI range 2-123 umol/L. Suri 2018 reports the column as 'Bilirubin' without the 'total' qualifier; in the related Mould-lab ADC papers the TBILI canonical convention is total bilirubin in umol/L (SI units), which matches the Table 1 reporting unit and range here.",
      source_name        = "BILIRUBIN"
    ),
    CREAT = list(
      description        = "Baseline serum creatinine concentration",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline. Power effect on MMAE clearance CL only (exponent -0.143). Reference 72.4 umol/L (Suri 2018 Table 1 overall-population mean). Pooled-population CREAT range 35-159 umol/L. Suri 2018 supplement Table S5 notes that creatinine concentration was used instead of creatinine clearance because creatinine clearance is collinear with body size, which is already in the model.",
      source_name        = "CREATININE"
    ),
    TUMTP_PCALCL = list(
      description        = "Primary cutaneous anaplastic large-cell lymphoma indicator (1 = pcALCL, 0 = other tumor types)",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = "Time-fixed. Categorical effect on ADC clearance only: cl_adc *= e_pcalcl_cl^TUMTP_PCALCL with e_pcalcl_cl = 0.728 (pcALCL patients have ~27% lower typical-value ADC CL than non-pcALCL patients). Reference category is 0 = non-pcALCL (the four other tumor types in the pooled cohort: Hodgkin lymphoma, systemic ALCL, mycosis fungoides, other CD30+ hematologic malignancies). 16 of 380 patients in the pooled cohort were pcALCL (4.2%).",
      source_name        = "PCALCL"
    ),
    ADA_POS = list(
      description        = "Anti-drug antibody positive in a 'newer-assay' study (1 = ever ADA-positive in NCT01990534 [Walewski 2016] or NCT01578499 [ALCANZA], 0 = otherwise)",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = "Time-varying once positive (Suri 2018 Methods: 'patients were treated as being positive at all times following the first time when ADA positivity was detected'). Multiplicative additive effect on ADC clearance: cl_adc *= (1 + e_adapos_cl * ADA_POS) with e_adapos_cl = 0.125 (12.5% higher CL when positive in a newer study). Mutually exclusive with ADA_POSOLD and ADA_MISSING. Reference category 0 = ADA-negative or not in a newer-assay study. The newer assay (sensitivity 23.573 ng/mL, drug tolerance 25 ug/mL) was used in NCT01990534 and the ALCANZA phase III trial.",
      source_name        = "ATAPOSNEW"
    ),
    ADA_POSOLD = list(
      description        = "Anti-drug antibody positive in an 'older-assay' study (1 = ever ADA-positive in NCT00430846, NCT00649584, NCT00848926, NCT00866047, 0 = otherwise)",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = "Time-varying once positive. Multiplicative additive effect on ADC clearance: cl_adc *= (1 + e_adaposold_cl * ADA_POSOLD) with e_adaposold_cl = 0.177 (17.7% higher CL when positive in an older study). Mutually exclusive with ADA_POS and ADA_MISSING. Reference category 0 = ADA-negative or not in an older-assay study. The older assay (sensitivity 4 ng/mL, drug tolerance 3,125 ng/mL) was used in NCT00430846, NCT00649584, NCT00848926, and NCT00866047. The 'older' vs 'newer' split is retained as a covariate because differing assay sensitivity / drug tolerance led to differing apparent ADA-on-CL effect sizes.",
      source_name        = "ATAPOSOLD"
    ),
    ADA_MISSING = list(
      description        = "Anti-drug antibody result missing (1 = ADA value missing, 0 = ADA value reported)",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = "Time-fixed. Multiplicative additive effect on ADC clearance: cl_adc *= (1 + e_adam_cl * ADA_MISSING) with e_adam_cl = 0.192 (19.2% higher typical-value CL when ADA result is missing). Mutually exclusive with ADA_POS and ADA_POSOLD. Reference category 0 = ADA result reported (positive or negative). 205 of 380 patients in the pooled cohort had missing ADA values (Suri 2018 Table S2).",
      source_name        = "ATAMISSING"
    )
  )

  population <- list(
    n_subjects     = 380L,
    n_studies      = 6L,
    age_range      = "12-87 years (overall median 37; ALCANZA subgroup median 61)",
    age_median     = "37 years",
    weight_range   = "39-168 kg",
    weight_median  = "76.46 kg (overall mean)",
    bsa_range      = "1.264-2.858 m^2",
    bsa_median     = "1.865 m^2 (overall mean)",
    sex_female_pct = 45,
    race_ethnicity = c(White = 83, Black = 6, Asian = 8, Other = 2),
    disease_state  = "CD30-positive malignancies pooled across six studies: Hodgkin lymphoma (HL, 246 patients = 64.7%), systemic ALCL (sALCL, 66 = 17.4%), mycosis fungoides (MF, 50 = 13.2%), primary cutaneous ALCL (pcALCL, 16 = 4.2%), and other hematologic malignancies (2 = 0.5%). The ALCANZA phase III subgroup (66 CTCL patients) is the focus of the exposure-response analyses, but the population PK model is fit to the full pooled 380-patient dataset.",
    dose_range     = "Phase I dose-ranging 0.1-1.8 mg/kg IV every 3 weeks (NCT00430846, NCT00649584); phase II/III 1.2 or 1.8 mg/kg IV every 3 weeks (NCT00848926, NCT00866047, NCT01990534, ALCANZA NCT01578499); 30-minute IV infusion. Licensed regimen is 1.8 mg/kg every 3 weeks with dose capping at 180 mg for patients >100 kg.",
    regions        = "Multinational (US, Europe, Asia).",
    study_phase    = "Phase I (NCT00430846 Younes 2010, NCT00649584 Fanale 2012); phase II (NCT00848926 Younes 2012 HL, NCT00866047 Pro 2012 sALCL, NCT01990534 Walewski 2016 HL); phase III (NCT01578499 ALCANZA Prince 2017 CTCL).",
    n_observations = "9,541 ADC + 9,669 MMAE concentration records (excluding ~3% BLOQ records) across 380 subjects; 22,660 total NONMEM records.",
    reference_subject = "BSA 1.865 m^2, ALB 36.81 g/L, TBILI 7.78 umol/L, CREAT 72.4 umol/L, TUMTP_PCALCL 0, ADA_POS 0, ADA_POSOLD 0, ADA_MISSING 0 (Suri 2018 Table 1 overall-population means; ADA-negative reference).",
    notes          = "Baseline characteristics from Suri 2018 Table 1 (overall column, n = 380). Suri 2018 supplement 1 (Methods) states continuous covariates were modeled as TVP = P_pop * (cov / cov_mean)^theta, normalized for the population mean — so the reference values used here are the Table 1 reported means rather than round numbers. The ALCANZA subgroup (n = 66 CTCL patients with MF or pcALCL) was older (median age 61 vs 37 overall), had higher albumin (mean 42.4 vs 36.8 g/L), and lower creatinine clearance (mean 105 vs 139 mL/min) than the overall population."
  )

  ini({
    # ADC structural parameters (Suri 2018 supplement Table S1; 3-compartment
    # linear with zero-order input and first-order elimination from central).
    lcl  <- log(0.0478); label("ADC clearance (CL, L/hr)")                                    # Suri 2018 Table S1 (sup 6): 0.0478 (2.7% RSE)
    lvc  <- log(3.5);    label("ADC central volume (V1, L)")                                  # Suri 2018 Table S1 (sup 6): 3.5 (1.2% RSE)
    lq   <- log(0.0673); label("ADC inter-compartmental clearance to peripheral 1 (Q2, L/hr)")# Suri 2018 Table S1 (sup 6): 0.0673 (3.1% RSE)
    lvp  <- log(3.67);   label("ADC peripheral volume 1 (V2, L)")                             # Suri 2018 Table S1 (sup 6): 3.67 (2.3% RSE)
    lq2  <- log(0.0125); label("ADC inter-compartmental clearance to peripheral 2 (Q3, L/hr)")# Suri 2018 Table S1 (sup 6): 0.0125 (3.3% RSE)
    lvp2 <- log(5.79);   label("ADC peripheral volume 2 (V3, L)")                             # Suri 2018 Table S1 (sup 6): 5.79 (1.3% RSE)

    # ADC covariate effects (Suri 2018 supplement Table S1; reference values
    # from supplement 1 statement "continuous covariates were normalized for
    # the population mean," with the means taken from Suri 2018 Table 1).
    e_bsa_vc       <- 1.27;   label("Power exponent of (BSA / 1.865) on ADC V1 (unitless)")              # Suri 2018 Table S1: 1.27 (4.9% RSE)
    e_bsa_cl       <- 0.457;  label("Power exponent of (BSA / 1.865) on ADC CL (unitless)")              # Suri 2018 Table S1: 0.457 (16.8% RSE)
    e_alb_cl       <- -0.496; label("Power exponent of (ALB / 36.81) on ADC CL (unitless)")              # Suri 2018 Table S1: -0.496 (3.6% RSE)
    e_pcalcl_cl    <- 0.728;  label("Power-form multiplier of pcALCL on ADC CL: cl *= e_pcalcl_cl^TUMTP_PCALCL") # Suri 2018 Table S1: 0.728 (8.9% RSE)
    e_adapos_cl    <- 0.125;  label("Multiplicative additive effect of ADA-positive (newer-assay study) on ADC CL: cl *= (1 + e_adapos_cl * ADA_POS)") # Suri 2018 Table S1: 0.125 (10.1% RSE)
    e_adaposold_cl <- 0.177;  label("Multiplicative additive effect of ADA-positive (older-assay study) on ADC CL: cl *= (1 + e_adaposold_cl * ADA_POSOLD)") # Suri 2018 Table S1: 0.177 (6.0% RSE)
    e_adam_cl      <- 0.192;  label("Multiplicative additive effect of ADA-missing on ADC CL: cl *= (1 + e_adam_cl * ADA_MISSING)") # Suri 2018 Table S1: 0.192 (9.4% RSE)

    # MMAE structural parameters (Suri 2018 supplement Table S3; 2-compartment
    # linear with an upstream Target binding pool and Lag compartment fed by
    # ADC, mirroring the Mould-lab ADC framework also used in Zhou 2025).
    lcl_mmae  <- log(0.577);   label("MMAE clearance (CLM, L/hr)")                                              # Suri 2018 Table S3 (sup 8): 0.577 (1.2% RSE)
    lvc_mmae  <- log(16.0);    label("MMAE central volume (VM, L)")                                             # Suri 2018 Table S3 (sup 8): 16.0 (1.4% RSE)
    lq_mmae   <- log(2.65);    label("MMAE inter-compartmental clearance (QM, L/hr)")                           # Suri 2018 Table S3 (sup 8): 2.65 (1.2% RSE)
    lvp_mmae  <- log(14.2);    label("MMAE peripheral volume (VMP, L)")                                         # Suri 2018 Table S3 (sup 8): 14.2 (1.1% RSE)
    lkd_mmae   <- log(0.00069); label("MMAE binding rate constant (Kd, 1/hr)")                                  # Suri 2018 Table S3 (sup 8): 0.00069 (1.6% RSE) — supplement table column header reads "Kd 1/h"; main paper page 994 mistakenly prints the unit as "L/h" (see vignette Errata section)
    lalfm_mmae <- log(2.64);    label("Decay rate of ADC->MMAE proteolytic-conversion fraction (ALFM, 1/hr)")   # Suri 2018 Table S3 (sup 8): 2.64 (1.0% RSE) — main paper page 994 mistakenly prints the unit as "L/h" (see vignette Errata section)
    lklag_mmae <- log(15.7);    label("Lag-compartment empty rate constant (Klag, 1/hr)")                       # Suri 2018 Table S3 (sup 8): 15.7 (1.0% RSE) — main paper page 994 mistakenly prints the unit as "L/h" (see vignette Errata section)
    # FM (fraction metabolized) is fixed to 1 in Suri 2018 Table S3 — encoded
    # as a literal constant in model() rather than an estimated parameter.

    # MMAE covariate effects (Suri 2018 supplement Table S3).
    e_bsa_cl_mmae    <- 2.81;   label("Power exponent of (BSA / 1.865) on MMAE CL (unitless)")             # Suri 2018 Table S3: 2.81 (6.4% RSE)
    e_bsa_vc_mmae    <- 0.89;   label("Power exponent of (BSA / 1.865) on MMAE VM (unitless)")             # Suri 2018 Table S3: 0.89 (9.5% RSE)
    e_alb_cl_mmae    <- 0.982;  label("Power exponent of (ALB / 36.81) on MMAE CL (unitless)")             # Suri 2018 Table S3: 0.982 (3.2% RSE)
    e_tbili_cl_mmae  <- -0.1;   label("Power exponent of (TBILI / 7.78) on MMAE CL (unitless)")            # Suri 2018 Table S3: -0.1 (7.5% RSE)
    e_creat_cl_mmae  <- -0.143; label("Power exponent of (CREAT / 72.4) on MMAE CL (unitless)")            # Suri 2018 Table S3: -0.143 (10.2% RSE)

    # IIV (log-normal). %CV from Suri 2018 Tables S1 and S3; converted via
    # omega^2 = log(CV^2 + 1). Only CL and V1 of each model carry IIV per
    # Tables S1 / S3 (other parameters: "-" indicates not estimated / no IIV).
    etalcl   ~ 0.1484         # Suri 2018 Table S1: ADC CL  40.0% CV; omega^2 = log(0.40^2 + 1)  = 0.1484
    etalvc   ~ 0.02197        # Suri 2018 Table S1: ADC V1  14.9% CV; omega^2 = log(0.149^2 + 1) = 0.02197

    etalcl_mmae  ~ 0.1664     # Suri 2018 Table S3: MMAE CL 42.5% CV; omega^2 = log(0.425^2 + 1) = 0.1664
    etalvc_mmae  ~ 0.3725     # Suri 2018 Table S3: MMAE V1 66.7% CV; omega^2 = log(0.667^2 + 1) = 0.3725

    # Residual variability. Suri 2018 supplement 1 specifies the LTBS (log
    # transform both sides) residual model — a homoscedastic additive error
    # on the log scale, which equals proportional error on the linear scale
    # in nlmixr2. %CV from Tables S1 and S3 maps directly to propSd.
    propSd        <- 0.291; label("Proportional residual error on ADC Cc (fraction)")                       # Suri 2018 Table S1: 29.1% CV (0.3% RSE)
    propSd_mmae     <- 0.423; label("Proportional residual error on MMAE Cc_mmae (fraction)")                 # Suri 2018 Table S3: 42.3% CV (0.3% RSE)
  })

  model({
    # 1. Derived covariate normalizations (Suri 2018 Table 1 overall-population
    # means used as reference values, per supplement 1's "normalized for the
    # population mean" convention).
    nbsa   <- BSA   / 1.865
    nalb   <- ALB   / 36.81
    ntbili <- TBILI / 7.78
    ncreat <- CREAT / 72.4

    # 2. Individual ADC parameters
    cl_adc <- exp(lcl + etalcl) *
              nbsa^e_bsa_cl *
              nalb^e_alb_cl *
              e_pcalcl_cl^TUMTP_PCALCL *
              (1 + e_adapos_cl    * ADA_POS) *
              (1 + e_adaposold_cl * ADA_POSOLD) *
              (1 + e_adam_cl      * ADA_MISSING)
    v1_adc <- exp(lvc + etalvc) * nbsa^e_bsa_vc
    q2_adc <- exp(lq)
    v2_adc <- exp(lvp)
    q3_adc <- exp(lq2)
    v3_adc <- exp(lvp2)

    # 2b. Individual MMAE parameters
    cl_mmae <- exp(lcl_mmae + etalcl_mmae) *
               nbsa^e_bsa_cl_mmae *
               nalb^e_alb_cl_mmae *
               ntbili^e_tbili_cl_mmae *
               ncreat^e_creat_cl_mmae
    vc_mmae <- exp(lvc_mmae + etalvc_mmae) * nbsa^e_bsa_vc_mmae
    qm      <- exp(lq_mmae)
    vp_mmae <- exp(lvp_mmae)
    kd      <- exp(lkd_mmae)
    alfm    <- exp(lalfm_mmae)
    klag    <- exp(lklag_mmae)

    # 3. Micro-constants
    k10 <- cl_adc  / v1_adc
    k12 <- q2_adc  / v1_adc
    k21 <- q2_adc  / v2_adc
    k13 <- q3_adc  / v1_adc
    k31 <- q3_adc  / v3_adc
    k40 <- cl_mmae / vc_mmae
    k47 <- qm      / vc_mmae
    k74 <- qm      / vp_mmae

    # 4. Time-varying ADC->MMAE proteolytic-conversion fraction. Suri 2018
    # describes this as "the fraction of MMAE formed directly from ADC was
    # assumed to decrease following ADC administration, relative to time
    # after dose," with ALFM the rate constant. FM is fixed at 1 in Table
    # S3, so the proteolytic flux at tad = 0 equals the ADC elimination flux
    # k10 * central and decays exponentially with rate ALFM. tad()
    # returns time after the most recent dose.
    fmt <- exp(-alfm * tad())

    # 5. ODE system (Suri 2018 Figure 1b; Mould-lab ADC framework as also
    # implemented in Zhou_2025_brentuximab.R). target(0) = 1 sets the
    # initial unitless target pool to one; binding ADC depletes it
    # irreversibly. The Lag compartment receives flux from both the
    # target-binding pathway (Kd*Target*ADC) and the proteolytic pathway
    # (FM*exp(-ALFM*tad)*K10*ADC), and empties to MMAE central with rate
    # Klag.
    target(0) <- 1
    d/dt(central)          <- -(k10 + k12 + k13) * central + k21 * peripheral1 + k31 * peripheral2
    d/dt(peripheral1)      <-  k12 * central - k21 * peripheral1
    d/dt(peripheral2)      <-  k13 * central - k31 * peripheral2
    d/dt(target)           <- -kd * target * central
    d/dt(lag)              <-  kd * target * central + fmt * k10 * central - klag * lag
    d/dt(central_mmae)     <-  klag * lag - (k40 + k47) * central_mmae + k74 * peripheral1_mmae
    d/dt(peripheral1_mmae) <-  k47 * central_mmae - k74 * peripheral1_mmae

    # 6. Observations and residual error. Concentrations are in umol/L
    # (= micromolar = uM). Convert to clinical units in post-processing:
    # Cc (umol/L) * MW_ADC_kDa (~153.4) = Cc in ug/mL (paper uses ug/mL
    # and ng/mL for ADC; ug/mL = mg/L). Cc_mmae (umol/L) * MW_MMAE_Da
    # (~718.0) = Cc_mmae in ng/mL (paper uses ng/mL and pg/mL for MMAE).
    # Doses provided in mg must be converted to umol before dosing this
    # model: amt_umol = dose_mg / MW_ADC_kDa.
    Cc      <- central      / v1_adc
    Cc_mmae <- central_mmae / vc_mmae

    Cc      ~ prop(propSd)
    Cc_mmae ~ prop(propSd_mmae)
  })
}
