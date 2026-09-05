Downes_2023_vancomycin_full <- function() {
  description <- "Two-compartment IV population PK model for vancomycin in critically ill children aged 1-17 years in the PICU (Downes 2023, full biomarker model). Clearance and intercompartmental clearance scale allometrically with body weight (exponent 0.75, reference 27 kg); clearance additionally scales as a power function of cystatin-C-based estimated GFR from the Hoek equation (exponent 0.85, reference 134 mL/min/1.73 m^2) and as an exponential function of the natural logarithm of creatinine-normalized urinary neutrophil gelatinase-associated lipocalin (base 0.94, no reference normalization); central and peripheral volumes scale linearly with body weight (exponent 1, reference 27 kg). This is the paper's full model and the only one of the three carrying a urinary kidney-injury biomarker; it did not outperform the two reduced sibling models for Bayesian AUC estimation, and its precision degraded markedly under limited sampling."
  reference <- "Downes KJ, Zuppa AF, Sharova A, Neely MN. Optimizing Vancomycin Therapy in Critically Ill Children: A Population Pharmacokinetics Study to Inform Vancomycin Area under the Curve Estimation Using Novel Biomarkers. Pharmaceutics. 2023;15(5):1336. doi:10.3390/pharmaceutics15051336"
  vignette <- "Downes_2023_vancomycin"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Vancomycin was given as an intermittent IV infusion, so
  # the dose enters `central` directly and there is no depot state. Both states
  # are left unverified: Downes 2023 Methods section 2.2 states that vancomycin
  # concentrations were measured by chemiluminescent microparticle immunoassay
  # in samples "collected via arterial catheter, peripheral venipuncture, or
  # venous catheter" but never names the matrix as serum or plasma, and the
  # paper never states what the peripheral distribution compartment represents.
  # "plasma" follows the repository default and is not a paper-sourced claim.
  compartmentData <- list(
    central     = list(analyte = "vancomycin", units = "mg", specimen = "plasma", verified = FALSE),
    peripheral1 = list(analyte = "vancomycin", units = "mg", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Downes 2023 Methods section 2.3: 'Clearance parameters were allometrically scaled for standardized weight to a power of 0.75 and volume parameters were scaled by standardized weight (power 1); weight was standardized by the median of the subjects' weights (27 kg).' Note that the 27 kg normalizing constant printed in the model equation does not equal the 25.9 kg model-training-group median in Table 1; the printed constant is the one this model uses. Table 1 training-group weight median 25.9 kg (IQR 13.9-41.8).",
      source_name        = "WT"
    ),
    CRCL = list(
      description        = "Estimated glomerular filtration rate from the cystatin-C-based Hoek equation, BSA-normalized",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Downes 2023 Methods section 2.3: 'GFR was estimated using ... pCysC alone based on the Hoek equation [32]' (Hoek FJ et al., Nephrol Dial Transplant 2003;18:2024-2031). Stored under canonical CRCL, which covers BSA-normalized GFR estimates regardless of the filtration marker; the assay form here is a cystatin-C-based estimate, chosen precisely BECAUSE creatinine is unreliable in critically ill children (Downes 2023 Introduction). Reference value 134 mL/min/1.73 m^2 is the normalizing constant printed in Table 2 and Table S4; it does not equal the Table 1 training-group median of 143 mL/min/1.73 m^2 at PK sampling (IQR 110-197), and the printed constant is the one this model uses. Half of the training group met the paper's augmented-renal-clearance definition of eGFR > 130 mL/min/1.73 m^2, so the cohort skews supranormal and the model carries little information about renal impairment.",
      source_name        = "HOEK"
    ),
    UNGALCR = list(
      description        = "Urinary neutrophil gelatinase-associated lipocalin normalized to urine creatinine",
      units              = "ng/mg creatinine",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Downes 2023 Methods section 2.3: 'Urinary biomarkers were normalized to urine creatinine (i.e., [biomarker]/[uCr]) to account for urine volume.' Enters the clearance equation as the paper's LNGAL term, defined in the Table 2 and Table S4 footnotes as 'the natural logarithm of the urinary NGAL concentration normalized to urinary creatinine (uNGAL/uCr)', so this column holds the RATIO on its natural scale and model() takes the logarithm. IMPORTANT: the term is NOT normalized to a reference value, so CL0 = 3.31 L/h is the clearance at UNGALCR = 1 ng/mg creatinine, which is far below any observed value; at the Table S3 training-group median of 60.3 ng/mg the NGAL factor is 0.94^log(60.3) = 0.776, a 22 percent clearance reduction. Table S3 training-group median 60.3 ng/mg creatinine (IQR 24.1-249.8); Discussion reports the training-group range as 1.4-2809 and the testing-group range as 7.8-10034, a spread the authors flag as possibly too imprecise for individual-level AUC estimation. Urine was collected on the evening and morning prior to PK sampling, with residual pre-vancomycin clinical samples used as baseline where available.",
      source_name        = "uNGAL/uCr"
    )
  )

  # Screened in the Downes 2023 covariate search (Table S2) but NOT retained in
  # the full model, so they are documentation only and are not referenced in
  # model(). Table S2 step 3 shows that urinary KIM-1, cystatin C and
  # osteopontin all improved the fit when log-transformed and creatinine-
  # normalized, but urinary NGAL gave the largest -2*LL and AIC reduction and
  # was the only urinary biomarker retained.
  covariatesDataExcluded <- list(
    AGE = list(
      description        = "Subject age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened on CL as a power function of (AGE/10) (Table S2 step 3, AIC change -1.9) and as a Hill function (AIC change -0.6), and on V1 as a power function (AIC change -2.7). Not retained.",
      source_name        = "AGE"
    ),
    SEXF = list(
      description        = "Female sex indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "male",
      notes              = "Screened on CL (Table S2 step 3, AIC change -4.0) and on V1 (AIC change -0.2), parameterized as TH2^FEM with FEM = 1 if female and 0 if male, the same orientation as canonical SEXF. Not retained.",
      source_name        = "FEM"
    ),
    CONMED_INOTROPE = list(
      description        = "Active vasopressor receipt indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "no active vasopressor receipt",
      notes              = "Screened on CL (Table S2 step 3, AIC change -5.9) and on V1 (AIC change -1.4), parameterized as TH2^VASO with VASO = 1 if active vasopressor receipt. Table 1 reports vasopressor receipt in 47 percent of the training group at PK sampling. Not retained. Stored under canonical CONMED_INOTROPE, which covers concomitant inotrope / vasoactive coadministration.",
      source_name        = "VASO"
    ),
    PIM3 = list(
      description        = "Pediatric Index of Mortality 3 score, expressed as a probability of death",
      units              = "(probability, 0-1)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened on CL and on V1 as (PIM3/2.82)^TH2 (Table S2 step 3, AIC changes -1.6 and -3.0). Table 1 training-group median 1.3 percent (IQR 0.5-4.3) at PK sampling. Not retained. No canonical register entry is proposed because the covariate is not used by any registered model; the name here is the paper's own and is documentation only.",
      source_name        = "PIM3"
    ),
    CREAT = list(
      description        = "Serum creatinine",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened on CL as (0.4/SCR)^TH1 (Table S2 step 2, AIC change -12.9; and again in step 3 with AIC change -3.6). Cystatin C outperformed serum creatinine by an AIC difference of -14.6, which is the paper's headline finding. Not retained in the full model; the sibling Downes_2023_vancomycin_schwartz model uses a creatinine-based eGFR instead.",
      source_name        = "SCR"
    ),
    CYSC = list(
      description        = "Plasma cystatin C",
      units              = "mg/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened on CL as (0.6/PCYSC)^TH1 (Table S2 step 2, AIC change -27.5), giving essentially the same fit as eGFR-Hoek (AIC change -27.2). Downes 2023 Results section 3.2 states the authors 'proceeded with further model training using the eGFRHoek model since clinical dosing guidance is typically based on a patient's eGFR rather than a direct biomarker result', so raw cystatin C is not in the final model even though it fit equally well. Table S3 training-group median 0.55 mg/L (IQR 0.4-0.7).",
      source_name        = "PCYSC"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 30L,
    n_studies        = 1L,
    n_sites          = 1L,
    n_concentrations = 150L,
    age_range        = "1-17 years (eligibility); Table 1 training-group median 9.8 years, IQR 3.8-11.2",
    age_median       = "9.8 years (IQR 3.8-11.2)",
    weight_range     = "IQR 13.9-41.8 kg (full range not reported)",
    weight_median    = "25.9 kg (IQR 13.9-41.8)",
    sex_female_pct   = 37,
    disease_state    = "Critically ill children in a single quaternary-care PICU receiving intermittent IV vancomycin for a suspected infection (a microbiological culture performed within 24 h of vancomycin initiation). Excluded: renal replacement therapy, plasmapheresis, extracorporeal membrane oxygenation, and age under 1 year (the authors excluded infants because cystatin C is affected by renal maturation over the first year of life, and state explicitly that the model cannot be applied to infants). PIM3 probability of death median 1.3 percent (IQR 0.5-4.3) at PK sampling; vasopressors in 47 percent at PK sampling.",
    dose_range       = "Clinician-chosen intermittent IV regimens; typical initial dosages 10-15 mg/kg/dose every 6-8 h. Table 1 training-group dose at PK sampling median 13.2 mg/kg (IQR 10.0-14.8)",
    regions          = "United States (Children's Hospital of Philadelphia, Philadelphia PA)",
    renal_function   = "eGFR-Hoek median 143 mL/min/1.73 m^2 (IQR 110-197) and eGFR-Schwartz median 164 (IQR 114-222) at PK sampling; serum creatinine median 0.30 mg/dL (IQR 0.20-0.48). Half of the training group had augmented renal clearance (eGFR-Hoek > 130 mL/min/1.73 m^2).",
    notes            = "Prospective observational study conducted August 2018 to July 2021. 50 evaluable subjects were split into a 30-subject model TRAINING group (this model) and a 20-subject model TESTING group used only for Bayesian AUC24 evaluation; the parameters here were estimated from the training group alone. The 30 training subjects contributed 150 vancomycin concentrations (14 clinical TDM samples and 136 research PK samples) spanning 3.9-67.8 ug/mL, with 2 values below the 3.0 ug/mL LLOQ coded as LLOQ/2. Model fit by NONPARAMETRIC adaptive grid (NPAG) in Pmetrics 1.9.7, not by a parametric method. The authors state explicitly (Discussion) that for a nonparametric fit 'the idea of typical parameter values and interindividual variability around them does not apply', so the Table 2 medians and CV percents encoded here are lognormal summaries of a joint density whose true shape is not reported -- see the vignette Assumptions section. Reported eta shrinkage is high across every parameter (49-61 percent, Table 2), consistent with only 150 concentrations across 30 subjects."
  )

  ini({
    # Structural parameters (Downes 2023 Table 2, replicated in Table S4 'Full
    # model' column). Values are weighted MEDIANS of the nonparametric joint
    # density. The reference subject weighs 27 kg, has an eGFR-Hoek of
    # 134 mL/min/1.73 m^2, and has UNGALCR = 1 ng/mg creatinine (see the
    # UNGALCR covariateData note -- the NGAL term carries no reference
    # normalization).
    lcl <- log(3.31); label("Clearance at WT=27 kg, CRCL=134 mL/min/1.73 m^2 and UNGALCR=1 ng/mg (L/h)")  # Table 2: CL0 3.31 (95th percentile 2.53-4.22)
    lvc <- log(3.50); label("Central volume at WT=27 kg (L)")                                             # Table 2: VC0 3.50 (95th percentile 2.72-7.09)
    lq  <- log(7.09); label("Intercompartmental clearance at WT=27 kg (L/h)")                             # Table 2: Q0 7.09 (95th percentile 4.76-7.97)
    lvp <- log(7.75); label("Peripheral volume at WT=27 kg (L)")                                          # Table 2: VP0 7.75 (95th percentile 6.63-13.80)

    # Allometric exponents. Downes 2023 Methods section 2.3 states the exponents
    # were imposed ("Clearance parameters were allometrically scaled ... to a
    # power of 0.75 and volume parameters were scaled by standardized weight
    # (power 1)"), and Table 2 prints CLWT / QWT / VC-WT / VP-WT with "-" in
    # both the 95th-percentile and CV% columns where every estimated parameter
    # carries a value, so they are fixed allometric constants.
    e_wt_cl_q  <- fixed(0.75); label("Allometric exponent on (WT/27) for CL and Q (unitless)")   # Table 2: CLWT = QWT = 0.75, no 95th percentile or CV%
    e_wt_vc_vp <- fixed(1);    label("Allometric exponent on (WT/27) for Vc and Vp (unitless)")  # Table 2: VC-WT = VP-WT = 1, no 95th percentile or CV%

    # Renal-function effect on clearance: (HOEK/134)^CLHOEK. Estimated.
    e_crcl_cl <- 0.85; label("Power exponent on (CRCL/134) for CL (unitless)")  # Table 2: CLHOEK 0.85 (95th percentile 0.22-0.90)

    # Urinary-NGAL effect on clearance. Note the FORM: the paper's Table 2
    # equation is (CLNGAL)^LNGAL, so CLNGAL is the BASE of an exponential whose
    # exponent is log(UNGALCR), not itself a power exponent. A value below 1
    # therefore means clearance FALLS as urinary NGAL rises, which is the
    # expected direction for a marker of kidney injury. Estimated.
    e_ungalcr_cl <- 0.94; label("Base of the 0.94^log(UNGALCR) factor on CL (unitless)")  # Table 2: CLNGAL 0.94 (95th percentile 0.86-1.00)

    # Interindividual variability. Downes 2023 Table 2 reports a CV% for every
    # ESTIMATED parameter of the nonparametric joint density -- including the
    # two covariate coefficients, because NPAG places every model parameter in
    # that density. Converted to lognormal variances with the standard identity
    # omega^2 = log(CV^2 + 1), which reproduces the reported CV% exactly and
    # preserves the reported value as the distribution median.
    #
    # The nonparametric joint density's CORRELATION structure is not reported,
    # so these etas are independent -- see the vignette Assumptions section.
    etalcl ~ 0.1449874  # log(0.395^2 + 1); Table 2: CL0    CV 39.5% [shrinkage 54.7%]
    etalvc ~ 0.2167745  # log(0.492^2 + 1); Table 2: VC0    CV 49.2% [shrinkage 59.5%]
    etalq  ~ 0.1009994  # log(0.326^2 + 1); Table 2: Q0     CV 32.6% [shrinkage 60.9%]
    etalvp ~ 0.1415864  # log(0.390^2 + 1); Table 2: VP0    CV 39.0% [shrinkage 49.1%]

    etae_crcl_cl    ~ 0.3324527  # log(0.628^2 + 1); Table 2: CLHOEK CV 62.8% [shrinkage 55.6%]
    etae_ungalcr_cl ~ 0.0118110  # log(0.109^2 + 1); Table 2: CLNGAL CV 10.9% [shrinkage 50.5%]

    # Residual variability. Downes 2023 Methods section 2.3: "The weighting
    # function on observations was 1/(gamma * SD^2), where SD ... was a combined
    # additive and multiplicative function of the assay imprecision as a
    # polynomial equation: SD = C0 + [C1 x observed concentration]; C0 had a
    # value of 1.5 (half the LLOQ) and C1 a value of 0.1, i.e., an assay with 10%
    # coefficient of variation (CV%). Gamma was initially set to 1 and fitted to
    # estimate the residual model error."
    #
    # C0 and C1 are stated assay constants, hence fixed(). The fitted gamma
    # multiplier is NOT reported anywhere in the paper or supplement, so it is
    # assumed to be 1 here; because gamma >= 1 in practice, the residual
    # magnitude encoded below is a LOWER BOUND on the model's true residual
    # error. Pmetrics forms the SD as a LINEAR sum of the additive and
    # proportional parts, which is nlmixr2's combined1() error structure, not
    # the default quadrature combination.
    addSd  <- fixed(1.5); label("Additive residual SD (ug/mL); assay C0, gamma assumed 1")           # Methods 2.3: C0 = 1.5 (half the 3.0 ug/mL LLOQ)
    propSd <- fixed(0.1); label("Proportional residual SD (fraction); assay C1, gamma assumed 1")    # Methods 2.3: C1 = 0.1 (10% assay CV)
  })
  model({
    # Individual PK parameters (Downes 2023 Table 2 footnote):
    #   CL = CL0 * (WT/27)^CLWT * (HOEK/134)^CLHOEK * (CLNGAL)^LNGAL
    #   V1 = VC0 * (WT/27)^VCWT
    #   Q  = Q0  * (WT/27)^QWT
    #   V2 = VP0 * (WT/27)^VPWT
    # where LNGAL is the natural logarithm of uNGAL/uCr.
    #
    # The two covariate coefficients carry their own lognormal IIV because the
    # NPAG joint density includes them; e * exp(eta) keeps the reported value as
    # the median and reproduces the reported CV% exactly.
    crcl_cl_i    <- e_crcl_cl    * exp(etae_crcl_cl)
    ungalcr_cl_i <- e_ungalcr_cl * exp(etae_ungalcr_cl)

    cl <- exp(lcl + etalcl) * (WT / 27)^e_wt_cl_q * (CRCL / 134)^crcl_cl_i * ungalcr_cl_i^log(UNGALCR)
    vc <- exp(lvc + etalvc) * (WT / 27)^e_wt_vc_vp
    q  <- exp(lq  + etalq)  * (WT / 27)^e_wt_cl_q
    vp <- exp(lvp + etalvp) * (WT / 27)^e_wt_vc_vp

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                  k12 * central - k21 * peripheral1

    # Dose in mg, volumes in L, so central/vc is mg/L == ug/mL.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd) + combined1()
  })
}
