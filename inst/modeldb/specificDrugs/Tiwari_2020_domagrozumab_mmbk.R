Tiwari_2020_domagrozumab_mmbk <- function() {
  description <- "Michaelis-Menten binding-kinetic (MM-BK) population PK/PD model for domagrozumab (anti-myostatin IgG1) and total myostatin in healthy adults (Tiwari 2020, the model selected by the authors): two-compartment IV/SC disposition with parallel linear and Michaelis-Menten elimination (Km decoupled from the binding constant), and a total-myostatin synthesis-degradation pool whose drug-bound fraction is cleared at the complex internalization rate through a quasi-steady-state binding constant Kss."
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
    # Structural PK, per kilogram (Tiwari 2020 Table 1, MM-BK column)
    lcl <- log(0.102); label("Linear clearance per kilogram (mL/h/kg)") # Table 1 MM-BK: CL = 0.102 mL/hour/kg (RSE 5.49%)
    lvc <- log(46); label("Central volume of distribution per kilogram (mL/kg)") # Table 1 MM-BK: VC = 46 mL/kg (RSE 11.6%)
    lq <- log(0.326); label("Intercompartmental clearance per kilogram (mL/h/kg)") # Table 1 MM-BK: Q = 0.326 mL/hour/kg (RSE 11.8%)
    lvp <- log(33); label("Peripheral volume of distribution per kilogram (mL/kg)") # Table 1 MM-BK: VP = 33 mL/kg (RSE 5.8%)
    lka <- log(0.009); label("First-order SC absorption rate constant (1/h)") # Table 1 MM-BK: ka = 0.009 1/hour (RSE 13.3%)
    lfdepot <- log(0.73); label("SC bioavailability (fraction)") # Table 1 MM-BK: Fbio = 0.73 (RSE 10.0%)

    # Saturable (Michaelis-Menten) elimination of free domagrozumab
    lvmax <- log(0.002); label("Maximum nonlinear elimination rate per kilogram (nmol/h/kg)") # Table 1 MM-BK: Vmax = 0.002 nmol/hour/kg (RSE 26%)
    lkm <- log(10.1); label("Free domagrozumab concentration at half-maximal nonlinear elimination (nmol/L)") # Table 1 MM-BK: KM = 10.1 nM (RSE 10.7%)

    # Total myostatin (Eq. 7)
    lkss <- log(4.39); label("Steady-state drug-myostatin binding constant (nmol/L)") # Table 1 MM-BK: KSS = 4.39 nM (RSE 3.99%)
    lkint <- log(0.009); label("First-order internalization rate of the drug-myostatin complex (1/h)") # Table 1 MM-BK: kint = 0.009 1/hour (RSE 7.14%)
    lkdeg <- log(0.046); label("First-order degradation rate of free myostatin (1/h)") # Table 1 MM-BK: kdeg = 0.046 1/hour (RSE 2.57%)
    lrbase <- log(0.149); label("Baseline total myostatin concentration M0 (nmol/L)") # Table 1 MM-BK: M0 = 0.149 nM (RSE 4.34%)

    # IIV, log-normal (Supplement 'Population Model'); Table 1 reports CV%.
    etalcl ~ log(1 + 0.213^2) # Table 1 MM-BK: IIV CL 21.3%
    etalvc ~ log(1 + 0.289^2) # Table 1 MM-BK: IIV VC 28.9%
    etalrbase ~ log(1 + 0.291^2) # Table 1 MM-BK: IIV M0 29.1%
    etalkint ~ log(1 + 0.273^2) # Table 1 MM-BK: IIV kint 27.3%
    etalvmax ~ log(1 + 0.621^2) # Table 1 MM-BK: IIV vmax 62.1%

    # Residual error: exponential (additive on log scale) for domagrozumab,
    # proportional for total myostatin (Supplement 'Population Model').
    expSd <- 0.184; label("Log-scale residual SD for free domagrozumab (log nmol/L)") # Table 1 MM-BK: RUV PK 18.4%
    propSd_Myo <- 0.22; label("Proportional residual error on total myostatin (fraction)") # Table 1 MM-BK: RUV PD 22.0%
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
    kss <- exp(lkss)
    kint <- exp(lkint + etalkint)
    kdeg <- exp(lkdeg)
    rbase <- exp(lrbase + etalrbase)
    ksyn <- kdeg * rbase # nmol/L/h; homeostasis at baseline

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Free domagrozumab in nmol/L: central (mg) / vc (mL) is g/L
    Cc_nM <- central / vc * 1e9 / mw_drug

    # Michaelis-Menten elimination as an amount rate: vmax (nmol/h) times
    # the saturation fraction, converted to mg/h
    mm_elim <- vmax * Cc_nM / (km + Cc_nM) * mw_drug * 1e-6

    # Eq. 6 (amount form) with Eqs. 2-3 for depot and peripheral
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - (kel + k12) * central + k21 * peripheral1 - mm_elim
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    f(depot) <- fdepot

    # Eq. 7: total myostatin; the drug-bound fraction Cc/(kss + Cc) is
    # cleared at kint instead of kdeg. Binding does not feed back on the PK.
    d/dt(total_target) <- ksyn - kdeg * total_target - (kint - kdeg) * total_target * Cc_nM / (kss + Cc_nM)
    total_target(0) <- rbase

    # Eq. 10: free myostatin and target coverage (% of myostatin bound)
    free_target <- total_target * kss / (kss + Cc_nM)
    coverage <- (1 - free_target / rbase) * 100

    Cc <- Cc_nM
    Myo <- total_target
    Cc ~ lnorm(expSd)
    Myo ~ prop(propSd_Myo)
  })
}
