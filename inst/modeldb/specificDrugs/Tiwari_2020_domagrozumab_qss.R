Tiwari_2020_domagrozumab_qss <- function() {
  description <- "Quasi-steady-state (QSS) target-mediated drug disposition population PK/PD model for domagrozumab (anti-myostatin IgG1) and total myostatin in healthy adults (Tiwari 2020): two-compartment IV/SC disposition of total domagrozumab with linear elimination of free drug, free drug solved from the QSS binding constant Kss, and a total-myostatin synthesis-degradation pool whose drug-bound complex is internalized at kint, which is also the only target-mediated elimination route of the drug."
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
    central = list(
      analyte = "domagrozumab (total: free + myostatin-bound)",
      units = "mg",
      specimen = "serum",
      verified = TRUE
    ),
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
    # Structural PK, per kilogram (Tiwari 2020 Table 1, QSS column)
    lcl <- log(0.108); label("Linear clearance of free drug per kilogram (mL/h/kg)") # Table 1 QSS: CL = 0.108 mL/hour/kg (RSE 4.56%)
    lvc <- log(48); label("Central volume of distribution per kilogram (mL/kg)") # Table 1 QSS: VC = 48 mL/kg (RSE 8.79%)
    lq <- log(0.398); label("Intercompartmental clearance per kilogram (mL/h/kg)") # Table 1 QSS: Q = 0.398 mL/hour/kg (RSE 44.7%)
    lvp <- log(18); label("Peripheral volume of distribution per kilogram (mL/kg)") # Table 1 QSS: VP = 18 mL/kg (RSE 13.5%)
    lka <- log(0.010); label("First-order SC absorption rate constant (1/h)") # Table 1 QSS: ka = 0.010 1/hour (RSE 11.7%)
    lfdepot <- log(0.62); label("SC bioavailability (fraction)") # Table 1 QSS: Fbio = 0.62 (RSE 8.0%)

    # Target binding and total myostatin (Eqs. 1, 4, 5)
    lkss <- log(1.99); label("Steady-state drug-myostatin binding constant (nmol/L)") # Table 1 QSS: KSS = 1.99 nM (RSE 2.72%)
    lkint <- log(0.013); label("First-order internalization rate of the drug-myostatin complex (1/h)") # Table 1 QSS: kint = 0.013 1/hour (RSE 8.25%)
    lkdeg <- log(0.063); label("First-order degradation rate of free myostatin (1/h)") # Table 1 QSS: kdeg = 0.063 1/hour (RSE 2.98%)
    lrbase <- log(0.141); label("Baseline total myostatin concentration M0 (nmol/L)") # Table 1 QSS: M0 = 0.141 nM (RSE 4.36%)

    # IIV, log-normal (Supplement 'Population Model'); Table 1 reports CV%.
    etalcl ~ log(1 + 0.196^2) # Table 1 QSS: IIV CL 19.6%
    etalvc ~ log(1 + 0.298^2) # Table 1 QSS: IIV VC 29.8%
    etalrbase ~ log(1 + 0.311^2) # Table 1 QSS: IIV M0 31.1%
    etalkint ~ log(1 + 0.312^2) # Table 1 QSS: IIV kint 31.2%

    # Residual error: exponential (additive on log scale) for domagrozumab,
    # proportional for total myostatin (Supplement 'Population Model').
    expSd <- 0.255; label("Log-scale residual SD for free domagrozumab (log nmol/L)") # Table 1 QSS: RUV PK 25.5%
    propSd_Myo <- 0.22; label("Proportional residual error on total myostatin (fraction)") # Table 1 QSS: RUV PD 22.0%
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
    kss <- exp(lkss)
    kint <- exp(lkint + etalkint)
    kdeg <- exp(lkdeg)
    rbase <- exp(lrbase + etalrbase)
    ksyn <- kdeg * rbase # nmol/L/h; homeostasis at baseline

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Total domagrozumab in the central compartment (nmol/L); central (mg)
    # holds free plus myostatin-bound drug
    Dtot_nM <- central / vc * 1e9 / mw_drug

    # Eq. 5: free domagrozumab from the QSS binding equilibrium
    bterm <- Dtot_nM - total_target - kss
    Cc_nM <- 0.5 * (bterm + sqrt(bterm^2 + 4 * kss * Dtot_nM))

    # Complex internalization (Eq. 1 last term) as an amount rate in mg/h:
    # nmol/L times vc (mL) times 1e-3 gives nmol; times mw_drug * 1e-6 gives mg
    complex_nM <- total_target * Cc_nM / (kss + Cc_nM)
    tmdd_elim <- kint * complex_nM * vc * mw_drug * 1e-9
    free_amt <- Cc_nM * vc * mw_drug * 1e-9 # mg of free drug in central

    # Eq. 1 (amount form): linear elimination and distribution act on free drug
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - (kel + k12) * free_amt + k21 * peripheral1 - tmdd_elim
    d/dt(peripheral1) <- k12 * free_amt - k21 * peripheral1
    f(depot) <- fdepot

    # Eq. 4: total myostatin
    d/dt(total_target) <- ksyn - kdeg * total_target - (kint - kdeg) * total_target * Cc_nM / (kss + Cc_nM)
    total_target(0) <- rbase

    # Free myostatin and target coverage, the Eq. 10 form applied to the QSS
    # model (same total-myostatin equation as the MM-BK model)
    free_target <- total_target * kss / (kss + Cc_nM)
    coverage <- (1 - free_target / rbase) * 100

    Cc <- Cc_nM
    Myo <- total_target
    Cc ~ lnorm(expSd)
    Myo ~ prop(propSd_Myo)
  })
}
