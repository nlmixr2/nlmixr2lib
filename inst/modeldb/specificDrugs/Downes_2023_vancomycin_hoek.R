Downes_2023_vancomycin_hoek <- function() {
  description <- "Two-compartment IV population PK model for vancomycin in critically ill children aged 1-17 years in the PICU (Downes 2023, reduced eGFR-Hoek model). Clearance and intercompartmental clearance scale allometrically with body weight (exponent 0.75, reference 27 kg) and clearance additionally scales as a power function of cystatin-C-based estimated GFR from the Hoek equation (exponent 1.00, reference 134 mL/min/1.73 m^2); central and peripheral volumes scale linearly with body weight (exponent 1, reference 27 kg). This is the reduced model the authors recommend for clinical Bayesian AUC estimation because cystatin C is available at many pediatric institutions, and it estimated AUC24 within 20 percent of the noncompartmental reference in 85 percent of subjects from a single optimally timed sample."
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
      notes              = "Downes 2023 Methods section 2.3: 'GFR was estimated using ... pCysC alone based on the Hoek equation [32]' (Hoek FJ et al., Nephrol Dial Transplant 2003;18:2024-2031). Stored under canonical CRCL, which covers BSA-normalized GFR estimates regardless of the filtration marker; the assay form here is a cystatin-C-based estimate, chosen precisely BECAUSE creatinine is unreliable in critically ill children (Downes 2023 Introduction). Reference value 134 mL/min/1.73 m^2 is the normalizing constant printed in Table S4 and in the Table S2 covariate-screening step; it does not equal the Table 1 training-group median of 143 mL/min/1.73 m^2 at PK sampling (IQR 110-197), and the printed constant is the one this model uses. Half of the training group met the paper's augmented-renal-clearance definition of eGFR > 130 mL/min/1.73 m^2, so the cohort skews supranormal and the model carries little information about renal impairment.",
      source_name        = "HOEK"
    )
  )

  # Screened in the Downes 2023 covariate search (Table S2) but NOT retained in
  # this reduced model, so they are documentation only and are not referenced in
  # model(). The full model (Downes_2023_vancomycin_full) retains urinary NGAL
  # in addition to eGFR-Hoek.
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
      notes              = "Screened on CL as (0.4/SCR)^TH1 (Table S2 step 2, AIC change -12.9; and again in step 3 with AIC change -3.6). Cystatin C outperformed serum creatinine by an AIC difference of -14.6, which is the paper's headline finding. Not retained.",
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
    notes            = "Prospective observational study conducted August 2018 to July 2021. 50 evaluable subjects were split into a 30-subject model TRAINING group (this model) and a 20-subject model TESTING group used only for Bayesian AUC24 evaluation; the parameters here were estimated from the training group alone. The 30 training subjects contributed 150 vancomycin concentrations (14 clinical TDM samples and 136 research PK samples) spanning 3.9-67.8 ug/mL, with 2 values below the 3.0 ug/mL LLOQ coded as LLOQ/2. Model fit by NONPARAMETRIC adaptive grid (NPAG) in Pmetrics 1.9.7, not by a parametric method, so the reported 'median' parameter values are weighted medians of the nonparametric joint density rather than typical values of an assumed distribution -- see the vignette Assumptions section for how that is represented here. Table S4 reports this reduced model as median and 95th-percentile range only, with no CV percent, so no interindividual variability could be encoded."
  )

  ini({
    # Structural parameters (Downes 2023 Table S4, eGFRHoek column). The
    # reference subject weighs 27 kg and has an eGFR-Hoek of
    # 134 mL/min/1.73 m^2.
    lcl <- log(2.61); label("Clearance at WT=27 kg and CRCL=134 mL/min/1.73 m^2 (L/h)")  # Table S4 eGFRHoek: CL0 2.61 (95th percentile 2.10-3.05)
    lvc <- log(5.09); label("Central volume at WT=27 kg (L)")                            # Table S4 eGFRHoek: VC0 5.09 (95th percentile 3.23-6.25)
    lq  <- log(6.18); label("Intercompartmental clearance at WT=27 kg (L/h)")            # Table S4 eGFRHoek: Q0 6.18 (95th percentile 4.47-7.76)
    lvp <- log(8.30); label("Peripheral volume at WT=27 kg (L)")                         # Table S4 eGFRHoek: VP0 8.30 (95th percentile 7.04-12.61)

    # Allometric exponents. Downes 2023 Methods section 2.3 states the exponents
    # were imposed ("Clearance parameters were allometrically scaled ... to a
    # power of 0.75 and volume parameters were scaled by standardized weight
    # (power 1)"), and Table 2 prints them without a CV percent or a
    # 95th-percentile range where every estimated parameter has both, so they
    # are fixed allometric constants rather than estimates.
    e_wt_cl_q  <- fixed(0.75); label("Allometric exponent on (WT/27) for CL and Q (unitless)")   # Table S4 equations: (WT/27)^0.75 on CL and on Q
    e_wt_vc_vp <- fixed(1);    label("Allometric exponent on (WT/27) for Vc and Vp (unitless)")  # Table S4 equations: (WT/27)^1 on V1 and on V2

    # Renal-function effect on clearance. Estimated, so not wrapped in fixed().
    e_crcl_cl <- 1.00; label("Power exponent on (CRCL/134) for CL (unitless)")  # Table S4 eGFRHoek: CLEGFR 1.00 (95th percentile 0.84-1.13)

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
    # Individual PK parameters (Downes 2023 Table S4, eGFRHoek model):
    #   CL = 2.61 * (WT/27)^0.75 * (eGFR/134)^CLEGFR
    #   V1 = 5.09 * (WT/27)^1
    #   Q  = 6.18 * (WT/27)^0.75
    #   V2 = 8.30 * (WT/27)^1
    # No eta terms: Table S4 reports this reduced model without CV percent.
    cl <- exp(lcl) * (WT / 27)^e_wt_cl_q * (CRCL / 134)^e_crcl_cl
    vc <- exp(lvc) * (WT / 27)^e_wt_vc_vp
    q  <- exp(lq)  * (WT / 27)^e_wt_cl_q
    vp <- exp(lvp) * (WT / 27)^e_wt_vc_vp

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
