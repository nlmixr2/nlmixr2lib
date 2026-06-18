Jiang_2008_gemcitabine <- function() {
  description <- paste(
    "Population PK model for intravenous gemcitabine and its primary",
    "inactive metabolite 2',2'-difluorodeoxyuridine (dFdU) in 94 adult",
    "patients with cancer pooled from three clinical studies (Jiang 2008",
    "Br J Clin Pharmacol 65(3):326-333). Both gemcitabine and dFdU are",
    "described by a two-compartment model with first-order elimination,",
    "joined by a first-order formation step. The fraction of gemcitabine",
    "converted to dFdU (F) is not identifiable, so dFdU parameters are",
    "apparent (CL/F, Q/F, Vc/F, Vp/F); the NONMEM ADVAN6 parent-metabolite",
    "encoding treats the total gemcitabine clearance as the apparent",
    "formation flux into the dFdU central compartment. Retained covariates",
    "after forward inclusion / backward elimination (Table 4): estimated",
    "creatinine clearance on apparent dFdU clearance via a linear-additive",
    "scaling CL_dFdU/F = 0.04 * (1 + 0.48 * CRCL/70); body surface area",
    "(power exponent 0.93, reference 1.73 m^2), oxaliplatin co-administration",
    "order (multiplicative factors 0.65 when gemcitabine is given first and",
    "0.54 when oxaliplatin is given first), and a non-small-cell lung cancer",
    "indicator (multiplicative factor 1.24) on apparent dFdU central volume.",
    "Residual error is proportional for both analytes."
  )
  reference <- paste(
    "Jiang X, Galettis P, Links M, Mitchell PL, McLachlan AJ. (2008).",
    "Population pharmacokinetics of gemcitabine and its metabolite in",
    "patients with cancer: effect of oxaliplatin and infusion rate.",
    "Br J Clin Pharmacol 65(3):326-333.",
    "doi:10.1111/j.1365-2125.2007.03040.x.",
    sep = " "
  )
  vignette <- "Jiang_2008_gemcitabine"
  units    <- list(time = "minute", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    CRCL = list(
      description        = "Estimated creatinine clearance (raw Cockcroft-Gault, NOT BSA-normalized).",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Linear-additive scaling on apparent dFdU clearance:",
        "CL_dFdU/F = 0.04 * (1 + 0.48 * CRCL/70) L/min (Jiang 2008 page 330",
        "covariate equation). The 70 mL/min denominator is a fixed",
        "normalisation constant taken from the source equation; the",
        "intercept at CRCL = 0 is 0.04 L/min. Median CRCL in the source",
        "cohort was 83 mL/min (range 37-150). Same orientation and",
        "assay-form precedent as Delattre 2010 amikacin (raw Cockcroft-",
        "Gault mL/min)."
      ),
      source_name        = "CLCR"
    ),
    BSA = list(
      description        = "Body surface area.",
      units              = "m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Power scaling on apparent dFdU central volume centred at 1.73",
        "m^2: V_C,dFdU/F = 46 * (BSA/1.73)^0.93 ... L (Jiang 2008 page",
        "330 covariate equation). Median BSA in the source cohort was 1.8",
        "m^2 (range 1.2-2.5). The paper text reports BSA was also",
        "significant on Q_dFdU/F, but the final covariate equation on",
        "page 330 retains BSA only on V_C,dFdU/F (see Table 4 forward-",
        "inclusion sequence)."
      ),
      source_name        = "BSA"
    ),
    GEMOX = list(
      description        = "Indicator for sequential gemcitabine-then-oxaliplatin combination dosing on the same study day.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "Gemcitabine alone (GEMOX = 0 AND OXGEM = 0)",
      notes              = paste(
        "1 = gemcitabine 30-min infusion followed by oxaliplatin 120-min",
        "infusion on the same day; 0 = otherwise. Multiplicative factor",
        "0.65 on V_C,dFdU/F when GEMOX = 1 (Jiang 2008 page 330 covariate",
        "equation: 0.65^GEMOX). Mutually exclusive with OXGEM: a patient",
        "is in exactly one of three groups -- gemcitabine alone (n=31,",
        "GEMOX=OXGEM=0), gem-then-ox (n=38, GEMOX=1), or ox-then-gem",
        "(n=25, OXGEM=1)."
      ),
      source_name        = "GEMOX"
    ),
    OXGEM = list(
      description        = "Indicator for sequential oxaliplatin-then-gemcitabine combination dosing on the same study day.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "Gemcitabine alone (GEMOX = 0 AND OXGEM = 0)",
      notes              = paste(
        "1 = oxaliplatin 120-min infusion followed by gemcitabine 30-min",
        "infusion on the same day; 0 = otherwise. Multiplicative factor",
        "0.54 on V_C,dFdU/F when OXGEM = 1 (Jiang 2008 page 330 covariate",
        "equation: 0.54^OXGEM). Mutually exclusive with GEMOX (see",
        "covariateData$GEMOX)."
      ),
      source_name        = "OXGEM"
    ),
    TUMTP_NSCLC = list(
      description        = "Non-small-cell lung cancer tumour-type indicator (1 = NSCLC, 0 = other).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "Other tumour types (pancreatic and diverse-tumour cohorts pooled).",
      notes              = paste(
        "Multiplicative factor 1.24 on V_C,dFdU/F when TUMTP_NSCLC = 1",
        "(Jiang 2008 page 330 covariate equation: 1.24^NSCLC). 47/94",
        "patients had NSCLC; the remaining 47 had pancreatic (n=14) or",
        "diverse (n=33) tumour types (Table 2). Same canonical as",
        "Ahamadi 2017 pembrolizumab and Aoyama 2012 sepantronium."
      ),
      source_name        = "NSCLC"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 94L,
    n_studies      = 3L,
    age_range      = "33-85 years",
    age_median     = "65 years",
    weight_range   = "37-120 kg",
    weight_median  = "73 kg",
    bsa_range      = "1.2-2.5 m^2",
    bsa_median     = "1.8 m^2",
    crcl_range     = "37-150 mL/min",
    crcl_median    = "83 mL/min",
    sex_female_pct = 27.7,
    disease_state  = paste(
      "Adult patients with cancer pooled from three clinical trials",
      "(Jiang 2008 Table 1). Study 1 (n=21) was an open-label Phase I",
      "dose-escalation in diverse tumour types receiving 30-min",
      "gemcitabine then 120-min oxaliplatin. Study 2 (n=43) was a",
      "randomised Phase II sequencing study in advanced NSCLC patients",
      "(gem-then-ox vs ox-then-gem). Study 3 (n=31) was a randomised",
      "crossover comparing 30-min vs 100-min gemcitabine infusion in",
      "diverse tumour types (no oxaliplatin). Tumour-type split: NSCLC",
      "47, pancreatic 14, diverse 33."
    ),
    dose_range     = paste(
      "Gemcitabine 750-1500 mg/m^2 as a 30-min IV infusion (or 1000",
      "mg/m^2 as a 100-min IV infusion in Study 3) with optional",
      "oxaliplatin 40-80 mg/m^2 as a 120-min IV infusion, on days 1 and",
      "8 of a 21-day cycle."
    ),
    regions        = "Australia (multi-centre Phase I/II oncology).",
    n_observations = paste(
      "652 gemcitabine and 1130 dFdU plasma concentrations across 122",
      "concentration-time profiles. Sampling at 0, 10, 25, 40 min and 1,",
      "1.5, 2, 2.5, 4, 7, 24 h and up to 192 h post gemcitabine dose.",
      "Limit of quantification 0.2 mg/L for both analytes (HPLC-UV at",
      "272 nm)."
    ),
    notes          = paste(
      "Demographics from Jiang 2008 Table 2. Software: NONMEM v5 level",
      "1.1 with ADVAN6 TRANS1 TOL5 and the FO estimation method.",
      "Covariate-development criterion DOFV >= 10.83 (P = 0.001, df = 1)",
      "for both forward inclusion and backward elimination, with a >20%",
      "clinical-significance threshold on parameter shifts. Age, body",
      "weight, sex and gemcitabine infusion rate were screened but not",
      "retained."
    )
  )

  ini({
    # Gemcitabine disposition: 2-compartment with first-order elimination
    # (Jiang 2008 Table 3). The full gemcitabine clearance flows into the
    # dFdU central compartment as an apparent formation flux because the
    # fraction-of-conversion F is not separately identifiable (NONMEM
    # ADVAN6 parent-metabolite parameterisation; see Hennig 2007
    # itraconazole / OH-itraconazole for the same idiom with f_m = 1).
    lcl <- log(2.7)  ; label("Gemcitabine systemic clearance CL_gem (L/min)")              # Jiang 2008 Table 3 CL_gem = 2.7 L/min
    lvc <- log(15)   ; label("Gemcitabine central volume V_C,gem (L)")                     # Jiang 2008 Table 3 V_C,gem = 15 L
    lq  <- log(0.7)  ; label("Gemcitabine inter-compartmental clearance Q_gem (L/min)")    # Jiang 2008 Table 3 Q_gem = 0.7 L/min
    lvp <- log(15)   ; label("Gemcitabine peripheral volume V_P,gem (L)")                  # Jiang 2008 Table 3 V_P,gem = 15 L

    # Apparent dFdU disposition (parameters scaled by 1/F): 2-compartment
    # with first-order elimination, formed at the rate (CL_gem/V_C,gem) *
    # gemcitabine central amount.
    lcl_dfdu <- log(0.04) ; label("Apparent dFdU clearance CL_dFdU/F at CRCL = 0 (L/min)")  # Jiang 2008 Table 3 CL_dFdU/F = 0.04 L/min; page 330 covariate equation CL_dFdU/F = 0.04 * (1 + 0.48 * CRCL/70)
    lvc_dfdu <- log(46)   ; label("Apparent dFdU central volume V_C,dFdU/F at BSA = 1.73 m^2, no oxaliplatin, no NSCLC (L)")  # Jiang 2008 Table 3 V_C,dFdU/F = 46 L; page 330 covariate equation V_C,dFdU/F = 46 * (BSA/1.73)^0.93 * 0.65^GEMOX * 0.54^OXGEM * 1.24^NSCLC
    lq_dfdu  <- log(0.2)  ; label("Apparent dFdU inter-compartmental clearance Q_dFdU/F (L/min)")  # Jiang 2008 Table 3 Q_dFdU/F = 0.2 L/min
    lvp_dfdu <- log(192)  ; label("Apparent dFdU peripheral volume V_P,dFdU/F (L)")        # Jiang 2008 Table 3 V_P,dFdU/F = 192 L

    # Covariate effects on the dFdU disposition (Jiang 2008 page 330).
    # The CRCL effect is linear-additive (1 + e * CRCL/70); BSA is a
    # power term centred at 1.73 m^2; the three binary indicators
    # (GEMOX, OXGEM, NSCLC) are multiplicative powers preserving the
    # source paper's symbolic form.
    e_crcl_cl_dfdu       <-  0.48  ; label("Linear slope of (CRCL/70) on apparent dFdU clearance (unitless)")  # Jiang 2008 page 330 covariate equation
    e_bsa_vc_dfdu        <-  0.93  ; label("Power exponent of (BSA/1.73) on apparent dFdU central volume (unitless)")  # Jiang 2008 page 330 covariate equation
    e_gemox_vc_dfdu      <-  0.65  ; label("Multiplicative factor on V_C,dFdU/F when gemcitabine is administered before oxaliplatin (unitless)")  # Jiang 2008 page 330 covariate equation: 0.65^GEMOX
    e_oxgem_vc_dfdu      <-  0.54  ; label("Multiplicative factor on V_C,dFdU/F when oxaliplatin is administered before gemcitabine (unitless)")  # Jiang 2008 page 330 covariate equation: 0.54^OXGEM
    e_tumtp_nsclc_vc_dfdu <- 1.24  ; label("Multiplicative factor on V_C,dFdU/F for NSCLC patients vs other tumour types (unitless)")  # Jiang 2008 page 330 covariate equation: 1.24^NSCLC

    # Inter-individual variability (BSV CV% from Jiang 2008 Table 3
    # converted to log-normal omega^2 via omega^2 = log(1 + CV^2)). The
    # paper also reports between-occasion variability (BOV CV%) on
    # CL_gem (12%), V_C,gem (93%), CL_dFdU/F (11%), and V_C,dFdU/F (10%);
    # nlmixr2 has no native occasion-level random-effect slot for a
    # static model definition, so BOV is documented as a deviation in
    # the validation vignette and only BSV is encoded here.
    etalcl       ~ 0.09175  # log(1 + 0.31^2); Jiang 2008 Table 3 BSV CL_gem = 31% CV
    etalvc       ~ 0.14165  # log(1 + 0.39^2); Jiang 2008 Table 3 BSV V_C,gem = 39% CV
    etalq        ~ 0.17699  # log(1 + 0.44^2); Jiang 2008 Table 3 BSV Q_gem = 44% CV
    etalvp       ~ 0.34336  # log(1 + 0.64^2); Jiang 2008 Table 3 BSV V_P,gem = 64% CV
    etalcl_dfdu  ~ 0.11553  # log(1 + 0.35^2); Jiang 2008 Table 3 BSV CL_dFdU/F = 35% CV
    etalvc_dfdu  ~ 0.02226  # log(1 + 0.15^2); Jiang 2008 Table 3 BSV V_C,dFdU/F = 15% CV (post-covariate residual; pre-oxaliplatin BSV was 20% per page 329)
    etalq_dfdu   ~ 0.08075  # log(1 + 0.29^2); Jiang 2008 Table 3 BSV Q_dFdU/F = 29% CV
    etalvp_dfdu  ~ 0.13489  # log(1 + 0.38^2); Jiang 2008 Table 3 BSV V_P,dFdU/F = 38% CV

    # Residual error: proportional for both gemcitabine and dFdU
    # (Jiang 2008 Methods: "Proportional statistical models were
    # employed to describe ... residual error for gemcitabine and dFdU").
    propSd      <- 0.40 ; label("Proportional residual error on gemcitabine (fraction)")  # Jiang 2008 Table 3 sigma_GEM = 40% CV
    propSd_dfdu <- 0.19 ; label("Proportional residual error on dFdU (fraction)")          # Jiang 2008 Table 3 sigma_dFdU = 19% CV
  })

  model({
    # Gemcitabine individual PK parameters.
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    q  <- exp(lq  + etalq)
    vp <- exp(lvp + etalvp)

    # Apparent dFdU individual PK parameters with the four covariate
    # effects on the typical values (Jiang 2008 page 330 covariate
    # equations). CRCL acts via a linear-additive scaling; BSA via a
    # power; GEMOX / OXGEM / TUMTP_NSCLC via multiplicative powers
    # preserving the source equations' symbolic form (a^X with X in
    # {0,1} is the multiplier when X = 1 and 1 when X = 0).
    cl_dfdu <- exp(lcl_dfdu + etalcl_dfdu) *
                 (1 + e_crcl_cl_dfdu * CRCL / 70)
    vc_dfdu <- exp(lvc_dfdu + etalvc_dfdu) *
                 (BSA / 1.73)^e_bsa_vc_dfdu *
                 e_gemox_vc_dfdu^GEMOX *
                 e_oxgem_vc_dfdu^OXGEM *
                 e_tumtp_nsclc_vc_dfdu^TUMTP_NSCLC
    q_dfdu  <- exp(lq_dfdu  + etalq_dfdu)
    vp_dfdu <- exp(lvp_dfdu + etalvp_dfdu)

    # Two-compartment gemcitabine disposition. All gemcitabine clearance
    # leaving the central compartment is treated as the apparent
    # formation flux into the dFdU central compartment (F-folded-in
    # NONMEM ADVAN6 parent-metabolite parameterisation; same idiom as
    # Hennig 2007 itraconazole with f_m fixed to 1). IV infusion is
    # delivered directly to central via the event-table rate/duration
    # column.
    d/dt(central)          <-  -(cl / vc) * central -
                                (q  / vc) * central +
                                (q  / vp) * peripheral1
    d/dt(peripheral1)      <-   (q  / vc) * central -
                                (q  / vp) * peripheral1

    # Two-compartment apparent dFdU disposition with first-order
    # formation from gemcitabine central at rate cl_gem / vc_gem.
    d/dt(central_dfdu)     <-   (cl / vc) * central -
                                (cl_dfdu / vc_dfdu) * central_dfdu -
                                (q_dfdu  / vc_dfdu) * central_dfdu +
                                (q_dfdu  / vp_dfdu) * peripheral1_dfdu
    d/dt(peripheral1_dfdu) <-   (q_dfdu  / vc_dfdu) * central_dfdu -
                                (q_dfdu  / vp_dfdu) * peripheral1_dfdu

    # Plasma concentration outputs.
    Cc      <- central      / vc
    Cc_dfdu <- central_dfdu / vc_dfdu

    Cc      ~ prop(propSd)
    Cc_dfdu ~ prop(propSd_dfdu)
  })
}
