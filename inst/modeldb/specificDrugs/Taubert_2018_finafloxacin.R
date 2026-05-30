Taubert_2018_finafloxacin <- function() {
  description <- paste(
    "Two-compartment population PK model for finafloxacin (a novel",
    "fluoroquinolone with enhanced antibacterial activity at acidic pH)",
    "with linear elimination, parallel first-order plus zero-order oral",
    "absorption (each with its own absorption lag time), an additive",
    "renal + non-renal clearance decomposition, and a cumulative-urinary",
    "excretion compartment. Built from pooled data of 266 subjects across",
    "three trials: 127 healthy volunteers (Trial I oral 25-1,000 mg/day;",
    "Trial II IV 200-1,000 mg/day) and 139 patients with complicated",
    "urinary tract infections (Trial III IV 800 mg/day, 60-min infusions).",
    "Covariates: body surface area on the central volume of distribution",
    "(power form, exponent 1.50, reference 1.829 m^2) and healthy /",
    "patient cohort status (DIS_HEALTHY) on both the renal and non-renal",
    "clearance arms. The paper-reported total apparent clearance (20.9 L/h",
    "healthy; -29% in patients) and population-specific fraction renally",
    "excreted (FER1 = 0.40 healthy, FER2 = 0.21 patient) are re-parameterised",
    "into the canonical lcl_renal + lcl_nonren additive decomposition;",
    "the typical values are anchored to DIS_HEALTHY = 0 (patient reference)",
    "per the inst/references/covariate-columns.md DIS_HEALTHY convention.",
    "The IIV translation between the paper and the re-parameterised forms",
    "is documented in the validation vignette Assumptions and deviations",
    "(Taubert 2018)."
  )
  reference <- paste(
    "Taubert M, Luckermann M, Vente A, Dalhoff A, Fuhr U.",
    "Population pharmacokinetics of finafloxacin in healthy volunteers",
    "and patients with complicated urinary tract infections.",
    "Antimicrob Agents Chemother. 2018;62(4):e02328-17.",
    "doi:10.1128/AAC.02328-17"
  )
  vignette <- "Taubert_2018_finafloxacin"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    BSA = list(
      description        = "Body surface area",
      units              = "m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Power scaling on Vc with reference BSA = 1.829 m^2 computed by",
        "the Mosteller formula (sqrt(W * H / 3600)) from the reference",
        "weight 70 kg and reference height 172 cm explicitly stated in",
        "Table 3 footnote a. The paper does not specify which BSA formula",
        "was used; Mosteller, DuBois, and Haycock all give ~1.83 m^2 at",
        "70 kg / 172 cm so the choice is numerically inconsequential for",
        "the reference value. Exponent on (BSA / 1.829) = 1.50 (Table 3).",
        "BSA was preferred over total body weight and lean body mass as",
        "the Vc scaling covariate (dOFV = -15.3 vs TBW)."
      ),
      source_name        = "BSA"
    ),
    DIS_HEALTHY = list(
      description        = "Healthy-participant cohort indicator: 1 = healthy adult volunteer (Trial I oral or Trial II IV), 0 = patient with complicated urinary tract infection (Trial III IV).",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = paste(
        "Time-fixed per subject. Paper Table 3 reports a fractional",
        "patient-status effect on total CL of -0.29 (CL_patient = 0.71 x",
        "CL_healthy) and two separate FER point estimates (FER1 = 0.40",
        "healthy, FER2 = 0.21 patient). The model file re-parameterises",
        "the paper's CL_total + FER decomposition into the canonical",
        "additive arms lcl_renal + lcl_nonren and re-anchors the typical",
        "to DIS_HEALTHY = 0 (patient reference) per the covariate-columns.md",
        "DIS_HEALTHY convention. The healthy multiplier on CL_renal is",
        "exp(e_dis_healthy_cl_renal) = 8.36 / 3.12 = 2.679 and on CL_nonren",
        "is 12.54 / 11.72 = 1.070; paper-reported typical clearances are",
        "reproduced exactly at both DIS_HEALTHY levels.",
        "Source column in the paper: PATIENT (1 = patient); DIS_HEALTHY = 1 - PATIENT."
      ),
      source_name        = "PATIENT"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 266L,
    n_studies      = 3L,
    age_range      = "19-90 years (Trial I 19-55 [median 42.5]; Trial II 19-63 [median 41.5]; Trial III 19-90 [median 61])",
    weight_range   = "50.9-140 kg (Trial I 61-97.7; Trial II 50.9-112.4; Trial III 64-140)",
    height_range   = "138-200 cm",
    sex_female_pct = 50.4,
    race_ethnicity = "Not explicitly reported in the article.",
    disease_state  = paste(
      "Trial I and II: healthy adult volunteers (single and multiple",
      "ascending doses, oral or IV). Trial III: hospitalised patients with",
      "complicated urinary tract infections or acute complicated /",
      "uncomplicated pyelonephritis (predominantly female, mostly older",
      "than 60 years)."
    ),
    dose_range     = paste(
      "Trial I: oral 25 to 1,000 mg/day (single and multiple ascending).",
      "Trial II: IV 200 to 1,000 mg/day. Trial III: IV 800 mg once daily",
      "for 5 or 10 days delivered as 60-min short-duration infusions."
    ),
    regions        = "Germany (Universities of Cologne, Kiel, and Giessen).",
    notes          = paste(
      "Pooled population PK analysis (NONMEM 7.3 FOCE-I; Perl-speaks-NONMEM",
      "bootstrap n = 1,000) of three phase I/II trials in healthy volunteers",
      "(Trials I oral and II IV) and cUTI patients (Trial III IV;",
      "ClinicalTrials.gov NCT01928433). 3,235 plasma samples (1,293 Trial I,",
      "1,531 Trial II, 414 Trial III) and 633 urine samples (496 healthy +",
      "137 patients) were modelled. Urine fit was good in healthy volunteers",
      "(R^2 = 0.92) but poor in patients (R^2 = 0.34); the IIV reported for",
      "FER2 (62%) reflects this. Plasma fits were good across trials",
      "(R^2 = 0.87-0.97). The paper notes a dose-dependent CL decrease at",
      "1,000 mg that was not implemented in the final model. The oral",
      "absorption parameters were fit sequentially after fixing the IV-fit",
      "CL, Vc, Q, Vp. IOV (inter-occasion variability) is reported in Table 3",
      "but not encoded in this model file -- a single random-effect set is",
      "used per parameter; see vignette Assumptions and deviations."
    )
  )

  ini({
    # ------------------------------------------------------------------------
    # STRUCTURAL CLEARANCE PARAMETERS  --  Taubert 2018 Table 3, re-parameterised.
    #
    # Paper Table 3 reports total apparent clearance CL = 20.9 L/h (95% CI
    # 18.5-23.4, IIV 54%, IOV 8%) with a fractional patient-status effect of
    # -0.29 (CL_patient = 0.71 x CL_healthy = 14.84 L/h) and two
    # population-specific fractions renally excreted: FER1 = 0.40 (healthy,
    # IIV 19%) and FER2 = 0.21 (patient, IIV 62%).
    #
    # This file re-parameterises the paper's structure into the canonical
    # additive decomposition lcl_renal + lcl_nonren so that the renal and
    # non-renal clearance components are first-class ini() parameters
    # (consistent with the cefepime example registered in parameter-names.md
    # under Multi-component clearance).
    #
    # Numerical equivalence (paper Table 3 typicals reproduced exactly):
    #   CL_renal_healthy  = 20.9       * 0.40  =  8.36  L/h
    #   CL_renal_patient  = 14.84      * 0.21  =  3.12  L/h
    #   CL_nonren_healthy = 20.9       * 0.60  = 12.54  L/h
    #   CL_nonren_patient = 14.84      * 0.79  = 11.72  L/h
    # The typicals are anchored to DIS_HEALTHY = 0 (patient) per the
    # covariate-columns.md DIS_HEALTHY convention; the healthy effect on
    # log scale is the log of the ratio healthy / patient.
    # ------------------------------------------------------------------------

    lcl_renal  <- log(3.12) ;    label("Renal CL at DIS_HEALTHY = 0 (patient typical, L/h)")
    # Paper Table 3: CL_total_patient (14.84 = 20.9 * 0.71) x FER2 (0.21) = 3.12 L/h
    lcl_nonren <- log(11.72) ;   label("Non-renal CL at DIS_HEALTHY = 0 (patient typical, L/h)")
    # Paper Table 3: CL_total_patient (14.84) x (1 - FER2) (0.79) = 11.72 L/h
    lvc        <- log(46.9) ;    label("Central volume of distribution at BSA reference 1.829 m^2 (L)")
    # Paper Table 3: Vc = 46.9 L (IIV 20%, IOV 8%)
    lq         <- log(2.80) ;    label("Intercompartmental clearance Q (L/h)")
    # Paper Table 3: Q = 2.80 L/h (IIV 57%, IOV 18%)
    lvp        <- log(43.1) ;    label("Peripheral volume of distribution Vp (L)")
    # Paper Table 3: Vp = 43.1 L (IIV 67%, IOV 22%)

    # Covariate effects
    e_bsa_vc                <- 1.50 ;          label("Power exponent of (BSA / 1.829) on Vc (unitless)")
    # Paper Table 3: BSA exponent on Vc = 1.50 (referenced to BSA at 70 kg / 172 cm = 1.829 m^2)
    e_dis_healthy_cl_renal  <- 0.985 ;         label("Log-scale effect of healthy status on CL_renal (unitless)")
    # = log(8.36 / 3.12) = log(2.679) = 0.985; re-derived from paper Table 3
    # CL_total_healthy (20.9) x FER1 (0.40) / (CL_total_patient (14.84) x FER2 (0.21))
    e_dis_healthy_cl_nonren <- 0.068 ;         label("Log-scale effect of healthy status on CL_nonren (unitless)")
    # = log(12.54 / 11.72) = log(1.070) = 0.068; re-derived from paper Table 3
    # CL_total_healthy (20.9) x (1 - FER1) (0.60) / (CL_total_patient (14.84) x (1 - FER2) (0.79))

    # ------------------------------------------------------------------------
    # ORAL ABSORPTION PARAMETERS  --  paper Table 3 oral data (sequential fit).
    #
    # Parallel first-order + zero-order absorption arms feed the central
    # compartment. The fraction of the oral dose entering the first-order
    # arm is f1st (logit-transformed because 0 < f1st < 1); the remainder
    # (1 - f1st) enters the zero-order arm with duration D0. Total oral
    # bioavailability F applies to the sum of the two arms; each arm has its
    # own absorption lag time.
    # ------------------------------------------------------------------------
    lfdepot   <- log(0.75) ;     label("Oral bioavailability F (unitless, log-scale)")
    # Paper Table 3 oral: F = 0.75 (IIV 33%, IOV 32%)
    logitf1st <- logit(0.77) ;   label("Logit fraction of oral dose entering first-order arm (f1st)")
    # Paper Table 3 oral: f1st = 0.77 (IIV 39%, IOV 36%)
    lka       <- log(6.61) ;     label("First-order oral absorption rate constant Ka (1/h)")
    # Paper Table 3 oral: Ka = 6.61 1/h (IIV 33%, IOV 168%)
    ltlag1    <- log(0.22) ;     label("First-order absorption lag time LAG1 (h)")
    # Paper Table 3 oral: LAG1 = 0.22 h (IIV 30%, IOV 28%)
    ld0       <- log(7.77) ;     label("Zero-order absorption duration D0 (h)")
    # Paper Table 3 oral: D0 = 7.77 h (IIV 31%, IOV 31%)
    ltlag2    <- log(0.54) ;     label("Zero-order absorption lag time LAG2 (h)")
    # Paper Table 3 oral: LAG2 = 0.54 h (IIV 33%, IOV 32%)

    # ------------------------------------------------------------------------
    # INTER-INDIVIDUAL VARIABILITY  --  log-normal omega^2 = log(CV^2 + 1).
    #
    # Paper Table 3 reports IIV on the paper's parameter set: CL_total (54%),
    # Vc (20%), Q (57%), Vp (67%), FER1 (19%), FER2 (62%), F (33%),
    # f1st (39%), Ka (33%), LAG1 (30%), D0 (31%), LAG2 (33%). Random effects
    # are treated as independent (the paper-reported bootstrap table does
    # not include off-diagonal covariances).
    #
    # The re-parameterised lcl_renal / lcl_nonren IIVs are approximations
    # to the paper's CL_total + FER stochastic model. Under independence of
    # the paper's CL_total and FER etas, log-additivity gives:
    #   var(log cl_renal)  ~ var(log CL_total) + var(log FER)
    # for healthy (FER1 IIV 19%): ~ log(1 + 0.54^2) + log(1 + 0.19^2) ~ 0.287
    # for patient (FER2 IIV 62%): ~ log(1 + 0.54^2) + log(1 + 0.62^2) ~ 0.571
    # The urine fit was poor in patients (R^2 = 0.34, paper Results); the
    # healthy-derived combined IIV (~57% CV, omega^2 = 0.287) is used as the
    # canonical value for etalcl_renal and the patient-population deviation
    # is documented in vignette Assumptions and deviations.
    #
    # For cl_nonren the paper's CL_total IIV (54%) dominates because
    # (1 - FER) varies modestly (60% healthy vs 79% patient); use it directly.
    # ------------------------------------------------------------------------
    etalcl_renal  ~ log(1 + 0.57 ^ 2)    # combined CL_total (54%) + FER1 (19%) IIV, healthy-derived
    etalcl_nonren ~ log(1 + 0.54 ^ 2)    # paper Table 3: IIV on total CL (cl_nonren-dominated)
    etalvc        ~ log(1 + 0.20 ^ 2)    # paper Table 3: Vc IIV = 20%
    etalq         ~ log(1 + 0.57 ^ 2)    # paper Table 3: Q IIV = 57%
    etalvp        ~ log(1 + 0.67 ^ 2)    # paper Table 3: Vp IIV = 67%

    # Oral absorption IIVs (sequential fit on Trial I oral data)
    etalfdepot    ~ log(1 + 0.33 ^ 2)    # paper Table 3 oral: F IIV = 33%
    etalogitf1st  ~ log(1 + 0.39 ^ 2)    # paper Table 3 oral: f1st IIV = 39%
    etalka        ~ log(1 + 0.33 ^ 2)    # paper Table 3 oral: Ka IIV = 33%
    etaltlag1     ~ log(1 + 0.30 ^ 2)    # paper Table 3 oral: LAG1 IIV = 30%
    etald0        ~ log(1 + 0.31 ^ 2)    # paper Table 3 oral: D0 IIV = 31%
    etaltlag2     ~ log(1 + 0.33 ^ 2)    # paper Table 3 oral: LAG2 IIV = 33%

    # ------------------------------------------------------------------------
    # RESIDUAL ERROR  --  paper Table 3 (combined additive + proportional on
    # the linear concentration scale). Paper-reported values are interpreted
    # as SDs (proportional values 0.24 plasma / 0.33 urine correspond to 24%
    # / 33% CV, consistent with popPK-reporting convention and the tight
    # bootstrap 95% CIs around 0.24 [0.22, 0.26]). See vignette Assumptions
    # and deviations for the SD-vs-variance interpretation note.
    #
    # The paper additionally reports oral-data residuals from the sequential
    # oral fit (additive 0.03 mg/L; proportional 0.14). The IV-data plasma
    # residuals are used here as the unified model default; users running
    # an oral-only simulation may switch propSd and addSd to the oral values.
    # ------------------------------------------------------------------------
    addSd            <- 0.001 ;   label("Additive plasma residual SD (mg/L)")
    # Paper Table 3 (IV data): Additive RUV plasma = 0.001 (CI 0.001, 0.002)
    propSd           <- 0.24 ;    label("Proportional plasma residual SD (fraction)")
    # Paper Table 3 (IV data): Proportional RUV plasma = 0.24 (CI 0.22, 0.26)
    addSd_urineAmt   <- 0.001 ;   label("Additive urine cumulative-amount residual SD (mg)")
    # Paper Table 3 (IV data): Additive RUV urine = 0.001 (CI 0.001, 0.154)
    propSd_urineAmt  <- 0.33 ;    label("Proportional urine cumulative-amount residual SD (fraction)")
    # Paper Table 3 (IV data): Proportional RUV urine = 0.33 (CI 0.29, 0.36)
  })

  model({
    # ---------------------------------------------------------------------
    # Reference value for BSA scaling on Vc.
    # ---------------------------------------------------------------------
    bsa_ref <- 1.829   # m^2; Mosteller(70 kg, 172 cm) = sqrt(70 * 172 / 3600)

    # ---------------------------------------------------------------------
    # Individual PK parameters with covariate effects.
    # ---------------------------------------------------------------------
    cl_renal  <- exp(lcl_renal  + etalcl_renal  + e_dis_healthy_cl_renal  * DIS_HEALTHY)
    cl_nonren <- exp(lcl_nonren + etalcl_nonren + e_dis_healthy_cl_nonren * DIS_HEALTHY)
    cl        <- cl_renal + cl_nonren
    vc        <- exp(lvc + etalvc) * (BSA / bsa_ref) ^ e_bsa_vc
    q         <- exp(lq  + etalq)
    vp        <- exp(lvp + etalvp)

    # Oral absorption (individual)
    fdepot <- exp(lfdepot + etalfdepot)
    f1st   <- expit(logitf1st + etalogitf1st)
    ka     <- exp(lka + etalka)
    tlag1  <- exp(ltlag1 + etaltlag1)
    d0     <- exp(ld0 + etald0)
    tlag2  <- exp(ltlag2 + etaltlag2)

    # Micro-constants
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ---------------------------------------------------------------------
    # ODE system.
    #
    # Two parallel oral absorption arms feed central:
    #   - depot  (first-order; rate Ka, lag LAG1, dose-fraction F * f1st)
    #   - depot2 (zero-order;  duration D0, lag LAG2, dose-fraction F * (1 - f1st))
    #
    # The zero-order arm uses rxode2's modeled-duration mechanism: an oral
    # dose record with rate = -2 to cmt = "depot2" is delivered as a
    # constant infusion at rate amt/D0 over D0 hours. depot2 has a fast
    # first-order drain (ktr_depot2 = 1000 1/h) so the infused mass passes
    # through to central with negligible residence time -- depot2's
    # quasi-steady-state amount is amt/(D0 * ktr_depot2), six orders of
    # magnitude smaller than the absorption window for D0 = 7.77 h, and
    # > 99.99% of the infused mass is delivered to central within the D0
    # absorption window. This idiom encodes a true zero-order input to
    # central without claiming f(central) (which would otherwise contaminate
    # IV doses that target central directly).
    #
    # IV doses target central directly; no f(central) is declared so the
    # IV dose is delivered intact.
    # ---------------------------------------------------------------------
    ktr_depot2 <- 1000.0   # 1/h; near-instantaneous depot2 -> central drain

    d/dt(depot)       <- -ka * depot
    d/dt(depot2)      <- -ktr_depot2 * depot2
    d/dt(central)     <-  ka * depot + ktr_depot2 * depot2 -
                          kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                  k12 * central - k21 * peripheral1

    # Urinary excretion compartment: cumulative renally excreted amount
    # (mass entering urine = (CL_renal / Vc) * central).
    d/dt(urine)       <- (cl_renal / vc) * central

    # Bioavailability and timing for the two oral absorption arms.
    f(depot)     <- fdepot * f1st
    alag(depot)  <- tlag1

    f(depot2)    <- fdepot * (1 - f1st)
    alag(depot2) <- tlag2
    dur(depot2)  <- d0

    # ---------------------------------------------------------------------
    # Observations.
    # Cc: plasma concentration in mg/L (dose in mg, vc in L -> mg/L).
    # urineAmt: cumulative urinary amount in mg (a derived output paired
    #   with propSd_urineAmt / addSd_urineAmt residual error).
    # ---------------------------------------------------------------------
    Cc       <- central / vc
    urineAmt <- urine

    Cc       ~ add(addSd) + prop(propSd)
    urineAmt ~ add(addSd_urineAmt) + prop(propSd_urineAmt)
  })
}
