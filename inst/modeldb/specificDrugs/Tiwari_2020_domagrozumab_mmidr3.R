Tiwari_2020_domagrozumab_mmidr3 <- function() {
  description <- "Michaelis-Menten indirect-response type 3 (MM-IDR3) population PK/PD model for domagrozumab (anti-myostatin IgG1) and total myostatin in healthy adults (Tiwari 2020): the two-compartment IV/SC PK with parallel linear and Michaelis-Menten elimination of the MM-BK model, and total myostatin described by an indirect response in which free domagrozumab stimulates myostatin synthesis through a sigmoid Smax function."
  reference <- "Tiwari A, Bhattacharya I, Chan PLS, Harnisch L. Comparing Model Performance in Characterizing the PK/PD of the Anti-Myostatin Antibody Domagrozumab. Clin Transl Sci. 2020;13(1):125-136. doi:10.1111/cts.12693"
  vignette <- "Tiwari_2020_domagrozumab"
  units <- list(time = "h", dosing = "mg", concentration = "nmol/L")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Not a covariate in the source model: every structural PK parameter in Tiwari 2020 Table 1 is reported per kilogram (CL and Q in mL/h/kg, VC and VP in mL/kg, Vmax in nmol/h/kg) because doses were given in mg/kg. WT multiplies those per-kg values in model() so that doses can be supplied in mg; with a dose of (mg/kg x WT) every concentration is independent of WT.",
      source_name = "WT"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "domagrozumab", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "domagrozumab", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "domagrozumab", units = "mg", specimen = "tissue", verified = TRUE),
    total_target = list(
      analyte = "myostatin (total: free + domagrozumab-bound)",
      units = "nmol/L",
      specimen = "serum",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = NA_integer_,
    n_studies = 1L,
    disease_state = "Healthy adult volunteers (first-in-human Phase 1 study of domagrozumab).",
    dose_range = "Single IV doses of 1, 3, 10, 20 and 40 mg/kg (2-h infusion); a single SC dose of 3 mg/kg; and 10 mg/kg IV (2-h infusion) every 2 weeks for 3 doses.",
    regions = "Not reported in Tiwari 2020.",
    notes = "Tiwari 2020 Methods ('Clinical study, PK, and myostatin data') describes the dose cohorts and defers subject numbers and demographics to Bhattacharya 2017 (Clin Pharmacol Drug Dev 7:484-497, doi:10.1002/cpdd.386), which is not on disk; n_subjects is therefore not filled. 21 domagrozumab and 1 myostatin observation were dropped from the analysis (Supplement, Figure S3 legend). LLOQ 0.2 nM (domagrozumab) and 0.04 nM (total myostatin)."
  )

  ini({
    # Structural PK, per kilogram (Tiwari 2020 Table 1, MM-IDR3 column)
    lcl <- log(0.100); label("Linear clearance per kilogram (mL/h/kg)") # Table 1 MM-IDR3: CL = 0.100 mL/hour/kg (RSE 6.47%)
    lvc <- log(46); label("Central volume of distribution per kilogram (mL/kg)") # Table 1 MM-IDR3: VC = 46 mL/kg (RSE 13.7%)
    lq <- log(0.318); label("Intercompartmental clearance per kilogram (mL/h/kg)") # Table 1 MM-IDR3: Q = 0.318 mL/hour/kg (RSE 13.1%)
    lvp <- log(34); label("Peripheral volume of distribution per kilogram (mL/kg)") # Table 1 MM-IDR3: VP = 34 mL/kg (RSE 6.9%)
    lka <- log(0.009); label("First-order SC absorption rate constant (1/h)") # Table 1 MM-IDR3: ka = 0.009 1/hour (RSE 14.4%)
    lfdepot <- log(0.72); label("SC bioavailability (fraction)") # Table 1 MM-IDR3: Fbio = 0.72 (RSE 10.4%)

    # Saturable (Michaelis-Menten) elimination of free domagrozumab
    lvmax <- log(0.002); label("Maximum nonlinear elimination rate per kilogram (nmol/h/kg)") # Table 1 MM-IDR3: Vmax = 0.002 nmol/hour/kg (RSE 26.7%)
    lkm <- log(10.8); label("Free domagrozumab concentration at half-maximal nonlinear elimination (nmol/L)") # Table 1 MM-IDR3: KM = 10.8 nM (RSE 10.4%)

    # Total myostatin, indirect response with stimulated synthesis (Eq. 9)
    lkdeg <- log(0.007); label("First-order degradation rate of total myostatin (1/h)") # Table 1 MM-IDR3: kdeg = 0.007 1/hour (RSE 3.31%)
    lrbase <- log(0.145); label("Baseline total myostatin concentration M0 (nmol/L)") # Table 1 MM-IDR3: M0 = 0.145 nM (RSE 4.94%)
    lsmax <- log(6.67); label("Maximum fold stimulation of myostatin synthesis (unitless)") # Table 1 MM-IDR3: Smax = 6.67 (RSE 7.5%)
    lsc50 <- log(124); label("Free domagrozumab concentration at half-maximal stimulation (nmol/L)") # Table 1 MM-IDR3: SC50 = 124 nM (RSE 33.7%)
    lhill <- log(0.646); label("Hill coefficient of the stimulation (unitless)") # Table 1 MM-IDR3: gamma = 0.646 (RSE 3.58%)

    # IIV, log-normal (Supplement 'Population Model'); Table 1 reports CV%.
    etalcl ~ log(1 + 0.216^2) # Table 1 MM-IDR3: IIV CL 21.6%
    etalvc ~ log(1 + 0.288^2) # Table 1 MM-IDR3: IIV VC 28.8%
    etalrbase ~ log(1 + 0.314^2) # Table 1 MM-IDR3: IIV M0 31.4%
    etalvmax ~ log(1 + 0.618^2) # Table 1 MM-IDR3: IIV vmax 61.8%
    etalsc50 ~ log(1 + 0.778^2) # Table 1 MM-IDR3: IIV SC50 77.8%
    etalsmax ~ log(1 + 0.227^2) # Table 1 MM-IDR3: IIV Smax 22.7%

    # Residual error: exponential (additive on log scale) for domagrozumab,
    # proportional for total myostatin (Supplement 'Population Model').
    expSd <- 0.183; label("Log-scale residual SD for free domagrozumab (log nmol/L)") # Table 1 MM-IDR3: RUV PK 18.3%
    propSd_Myo <- 0.205; label("Proportional residual error on total myostatin (fraction)") # Table 1 MM-IDR3: RUV PD 20.5%
  })
  model({
    # Domagrozumab molecular weight, needed to convert mg doses to the
    # paper's nmol scale. Tiwari 2020 does not print it; 150,000 g/mol is the
    # value used for domagrozumab in Nguyen 2020 (Table 2, MW_Drug).
    mw_drug <- 150000

    # Individual parameters; per-kg values times WT give per-subject values
    cl <- exp(lcl + etalcl) * WT # mL/h
    vc <- exp(lvc + etalvc) * WT # mL
    q <- exp(lq) * WT # mL/h
    vp <- exp(lvp) * WT # mL
    ka <- exp(lka)
    fdepot <- exp(lfdepot)
    vmax <- exp(lvmax + etalvmax) * WT # nmol/h
    km <- exp(lkm)
    kdeg <- exp(lkdeg)
    rbase <- exp(lrbase + etalrbase)
    smax <- exp(lsmax + etalsmax)
    sc50 <- exp(lsc50 + etalsc50)
    hill <- exp(lhill)
    ksyn <- kdeg * rbase # nmol/L/h; homeostasis at baseline

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Free domagrozumab in nmol/L: central (mg) / vc (mL) is g/L
    Cc_nM <- central / vc * 1e9 / mw_drug

    # Michaelis-Menten elimination as an amount rate: vmax (nmol/h) times
    # the saturation fraction, converted to mg/h
    mm_elim <- vmax * Cc_nM / (km + Cc_nM) * mw_drug * 1e-6

    # PK module of the MM-BK model (Eq. 6, amount form; Eqs. 2-3)
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - (kel + k12) * central + k21 * peripheral1 - mm_elim
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    f(depot) <- fdepot

    # Eq. 9: domagrozumab stimulates total-myostatin synthesis. The
    # concentration is floored at 0 because a fractional power of a
    # solver-noise negative value is NaN.
    cpos <- max(Cc_nM, 0)
    stim <- smax * cpos^hill / (sc50^hill + cpos^hill)
    d/dt(total_target) <- ksyn * (1 + stim) - kdeg * total_target
    total_target(0) <- rbase

    Cc <- Cc_nM
    Myo <- total_target
    Cc ~ lnorm(expSd)
    Myo ~ prop(propSd_Myo)
  })
}
