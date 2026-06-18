McLaughlin_2024_tld_1 <- function() {
  description <- paste(
    "Joint parent (entrapped) + free + metabolite (doxorubicinol)",
    "population PK model for TLD-1, a novel small-diameter pegylated",
    "liposomal doxorubicin, in 30 adults with advanced solid tumours",
    "(McLaughlin 2024 phase I dose-escalation, SAKK 65/16 / NCT03387917).",
    "Structure: one-compartment liposomal-entrapped reservoir (V1) with",
    "linear release into the free-doxorubicin central compartment",
    "(release rate krel = CL1/V1); free doxorubicin disposition is two-",
    "compartment (Vc, Vp, Q) with linear metabolism-to-doxorubicinol",
    "clearance (CL); doxorubicinol is one-compartment with linear",
    "elimination (Vc_doxol, CL_doxol). Body surface area (BSA, reference",
    "1.75 m^2) is the only retained covariate, entering as a power",
    "model on the free-doxorubicin central (V2 exponent 4.47) and",
    "peripheral (V3 exponent 11.5) volumes. Inter-individual variability",
    "is fitted on CL1 (release), V1 (shared-eta scale 0.643 of",
    "ome_CL1), CL2 (free->doxol) and CL4 (doxol elimination);",
    "inter-occasion variability on CL1, CL2, V1, V2 from Table 2 is",
    "documented but not encoded structurally here (nlmixr2lib has no",
    "canonical occasion-column convention; see Hempel 2003 / Hong 2006",
    "precedents). Residual error is log-transformed-both-sides additive",
    "on the log scale -- equivalent to proportional in nlmixr2's linear",
    "space and encoded here as separate propSd per analyte. Distinct",
    "from Hempel 2003 (paediatric liposomal daunorubicin, total drug",
    "only) and Varatharajan 2016 (free daunorubicin + daunorubicinol",
    "in adult AML, no liposomal reservoir)."
  )
  reference <- paste(
    "Mc Laughlin AM, Hess D, Michelet R, Colombo I, Haefliger S, Bastian S,",
    "Rabaglio M, Schwitter M, Fischer S, Eckhardt K, Hayoz S, Kopp C, Klose",
    "M, Sessa C, Stathis A, Halbherr S, Huisinga W, Joerger M, Kloft C.",
    "Population pharmacokinetics of TLD-1, a novel liposomal doxorubicin,",
    "in a phase I trial. Cancer Chemother Pharmacol. 2024;94(3):349-360.",
    "doi:10.1007/s00280-024-04679-z.",
    sep = " "
  )
  vignette <- "McLaughlin_2024_tld_1"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  paper_specific_compartments <- c("entrapped")

  covariateData <- list(
    BSA = list(
      description        = "Body surface area at baseline. Enters as a power covariate on the free-doxorubicin central (V2) and peripheral (V3) volumes only; the paper screened BSA on V1 (entrapped) but the effect was not retained (Mc Laughlin 2024 Results / Covariate model paragraph).",
      units              = "m^2",
      type               = "continuous",
      reference_category = "n/a -- used with power scaling (BSA / 1.75)^exponent. Reference BSA is the cohort median 1.75 m^2 (Mc Laughlin 2024 Results / Clinical data paragraph).",
      notes              = paste(
        "Cohort median BSA 1.75 m^2 (range 1.44-2.44 m^2); 30 adults with",
        "advanced solid tumours, 80% female. The power exponents on V2",
        "(4.47) and V3 (11.5) are unusually large -- the paper notes the",
        "small central volume of free doxorubicin (V2 ~ 0.5 L) makes the",
        "BSA scaling steep but the effect was retained because it",
        "significantly improved model fit (Mc Laughlin 2024 Results /",
        "Covariate model paragraph)."
      ),
      source_name        = "BSA"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Subject age at baseline (screened in stepwise selection but not retained).",
      units       = "years",
      type        = "continuous",
      notes       = "Mc Laughlin 2024 Methods / Covariate submodel: 'Potential patient characteristics to be implemented as covariates in the model were pre-selected based on plausibility, previous reports, and availability in the dataset.' Results: 'No other covariates were identified' (only BSA retained).",
      source_name = "AGE"
    ),
    WT = list(
      description = "Body weight at baseline (screened, not retained; replacing BSA with WT did not improve fit).",
      units       = "kg",
      type        = "continuous",
      notes       = "Mc Laughlin 2024 Results / Covariate model: 'Replacing BSA with other body size descriptors, including body weight or lean body weight, did not improve model fit.'",
      source_name = "WT"
    ),
    HT = list(
      description = "Body height at baseline (screened, not retained).",
      units       = "cm",
      type        = "continuous",
      notes       = "Mc Laughlin 2024 Methods / Analysis dataset generation: 'Patient characteristics age, body weight, body height, and BSA ... were included in the dataset and available for testing.'",
      source_name = "HT"
    ),
    BMI = list(
      description = "Body mass index calculated from WT and HT (screened, not retained).",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = "Mc Laughlin 2024 Methods: 'additional body size descriptors lean body weight and body mass index (BMI) were calculated and available for testing as potential covariates.'",
      source_name = "BMI"
    ),
    LBW = list(
      description = "Lean body weight by Janmahasatian formula (screened, not retained).",
      units       = "kg",
      type        = "continuous",
      notes       = "Mc Laughlin 2024 Methods: lean body weight calculated per Janmahasatian 2005 (reference 24); 'No other covariates were identified.'",
      source_name = "LBW"
    ),
    SCR = list(
      description = "Serum creatinine (screened, not retained).",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Mc Laughlin 2024 Methods: serum creatinine and CKD-EPI eGFR were available for covariate testing; not retained.",
      source_name = "SCR"
    ),
    EGFR = list(
      description = "Estimated glomerular filtration rate by CKD-EPI (screened, not retained).",
      units       = "mL/min/1.73 m^2",
      type        = "continuous",
      notes       = "Mc Laughlin 2024 Methods reference 23 (Levey 2009 CKD-EPI).",
      source_name = "EGFR"
    ),
    ALT = list(
      description = "Alanine aminotransferase (screened, not retained).",
      units       = "U/L",
      type        = "continuous",
      notes       = "Mc Laughlin 2024 Methods / Analysis dataset generation.",
      source_name = "ALT"
    ),
    AST = list(
      description = "Aspartate aminotransferase (screened, not retained).",
      units       = "U/L",
      type        = "continuous",
      notes       = "Mc Laughlin 2024 Methods / Analysis dataset generation.",
      source_name = "AST"
    ),
    ALP = list(
      description = "Alkaline phosphatase (screened, not retained).",
      units       = "U/L",
      type        = "continuous",
      notes       = "Mc Laughlin 2024 Methods / Analysis dataset generation.",
      source_name = "ALP"
    ),
    BILI = list(
      description = "Total bilirubin (screened, not retained).",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Mc Laughlin 2024 Methods / Analysis dataset generation.",
      source_name = "BILI"
    ),
    SEXF = list(
      description = "Female sex indicator (cohort 80% female; sex not retained as a covariate).",
      units       = "binary",
      type        = "binary",
      notes       = "Mc Laughlin 2024 Methods / Covariate submodel describes 'Categorical covariates (e.g., sex) were implemented using fractional change models'; Results: 'No other covariates were identified.'",
      source_name = "SEXF"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 30L,
    n_studies       = 1L,
    age_range       = "38-83 years",
    age_median      = "67.5 years",
    bsa_range       = "1.44-2.44 m^2",
    bsa_median      = "1.75 m^2",
    bmi_range       = "16.5-42.2 kg/m^2",
    bmi_median      = "24.7 kg/m^2",
    lbw_range       = "38.7-77.6 kg",
    lbw_median      = "48.5 kg",
    sex_female_pct  = 80,
    disease_state   = paste(
      "Advanced solid tumours; breast cancer 43.3%, ovarian 20.0%,",
      "gastrointestinal 3.3%, other 33.3% of the cohort."
    ),
    dose_range      = paste(
      "TLD-1 IV infusion every 21 days, dose levels 1-7 = 10, 16, 23,",
      "30, 35, 40, 45 mg/m^2 of body surface area. Dose levels 1-6",
      "infused over 60 min; dose level 7 over 90 min. Up to 6 cycles",
      "(anthracycline-pretreated) or 9 cycles (anthracycline-naive).",
      "Dose level 6 (40 mg/m^2) was the recommended dose for the",
      "expansion cohort. TLD-1 = pegylated liposomal doxorubicin with",
      "average liposome diameter ~36 nm (Caelyx ~70 nm); 99% of",
      "doxorubicin entrapped at time of infusion per the manufacturer."
    ),
    regions         = "Switzerland (4 phase I centres, Swiss Group for Clinical Cancer Research SAKK 65/16)",
    clinical_trial  = "NCT03387917",
    n_concentrations = "1870 (n = 624 total doxorubicin, n = 623 free doxorubicin, n = 623 doxorubicinol)",
    sampling_window = paste(
      "Cycles 1 and 2; pre-dose, mid-infusion (0.5 h DL 1-6 / 0.75 h",
      "DL 7), end of infusion (1 h DL 1-6 / 1.5 h DL 7), 0.5, 1, 3, 5,",
      "7 h after end of infusion, 24, 48 (cycle 1 only), 168 (day 8),",
      "336 (day 15) h after infusion."
    ),
    assay           = paste(
      "LC-MS/MS by Swiss BioQuant AG, Reinach. Total doxorubicin",
      "(doxorubicin_entrapped + doxorubicin_free) calibration 20.0-20000",
      "ng/mL (inter-batch precision 3.2-6.3%, accuracy 100.0-104.3%);",
      "free doxorubicin 2.00-2000 ng/mL (precision 4.8-10.9%, accuracy",
      "99.0-101.3%); doxorubicinol 0.500-500 ng/mL (precision 5.3-12.1%,",
      "accuracy 98.0-99.0%). Entrapped concentration computed by",
      "subtracting free from total; entrapped fraction at start of",
      "infusion fixed at 100% based on in-house manufacturer data showing",
      ">99% entrapment."
    ),
    notes           = paste(
      "Baseline demographics from Mc Laughlin 2024 Results /",
      "Clinical data paragraph and Table 2 footnote b (reference BSA",
      "1.75 m^2). NONMEM 7.4 with FOCE+I; SIR for parameter precision.",
      "12 patients in dose-escalation cohorts (one to two doses each",
      "across DL 1-7); 9 additional patients in DL 7 expansion (45",
      "mg/m^2 was the tentative MTD); 9 additional patients in DL 6",
      "(40 mg/m^2) expansion after late cumulative toxicities at DL 7."
    )
  )

  ini({
    # ---------------------------------------------------------------
    # Structural parameters. Reference values are the final NONMEM
    # estimates from Mc Laughlin 2024 Table 2 ('Parameter estimates for
    # the final joint parent-metabolite population pharmacokinetic
    # model of liposomal doxorubicin, free doxorubicin and
    # doxorubicinol'). Volumes in L, clearances in L/h, rate constants
    # in 1/h; doses in mg, concentrations in mg/L.
    #
    # Liposomal reservoir (entrapped doxorubicin): one-compartment with
    # release-rate kinetics. V1 (apparent volume of the entrapped state,
    # = plasma volume per the paper Methods Structural submodel) and
    # CL1 (linear release clearance of the entrapped state) are both
    # reported in Table 2; the model below carries V1 as `vd` (paper-
    # named volume of distribution for the paper-mechanistic state) and
    # the release rate constant krel = CL1 / V1 (= 0.0271 / 3.39 =
    # 0.007994 1/h) as a derived rate. This preserves the published
    # leakage half-life t1/2,leak = log(2) / krel = 86.7 h (Mc Laughlin
    # 2024 Discussion paragraph 4) and avoids introducing a non-
    # canonical 'lcl_release' clearance term.
    # ---------------------------------------------------------------
    lvd      <- log(3.39)              ; label("Apparent volume of distribution of entrapped doxorubicin (L)")                               # Table 2 row V1, RSE 7%
    lkrel    <- log(0.0271 / 3.39)     ; label("First-order release rate constant of doxorubicin from the liposomal reservoir (1/h)")        # Table 2 row CL1 / V1 derivation; krel = CL1/V1 = 0.0271/3.39 = 0.007994 1/h (leakage t1/2 = 86.7 h)

    # Free doxorubicin (parent of the canonical central + peripheral
    # disposition in this file): two-compartment with linear metabolism
    # to doxorubicinol. V2 and V3 are reference values at BSA = 1.75
    # m^2; covariate scaling applied in model() below.
    lvc      <- log(0.531)             ; label("Central volume of distribution of free doxorubicin V2 at BSA 1.75 m^2 (L)")                  # Table 2 row V2, RSE 16%
    lvp      <- log(61.3)              ; label("Peripheral volume of distribution of free doxorubicin V3 at BSA 1.75 m^2 (L)")               # Table 2 row V3, RSE 19%
    lq       <- log(0.136)             ; label("Inter-compartmental clearance of free doxorubicin Q (L/h)")                                  # Table 2 row Q, RSE 18%
    lcl      <- log(0.450)             ; label("Linear clearance of free doxorubicin to doxorubicinol CL2 (L/h)")                            # Table 2 row CL2, RSE 11%

    # Doxorubicinol metabolite (canonical metabolite suffix _doxol):
    # one-compartment with linear elimination.
    lvc_doxol <- log(8152)             ; label("Volume of distribution of doxorubicinol V4 (L)")                                              # Table 2 row V4, RSE 12%
    lcl_doxol <- log(74.6)             ; label("Linear clearance of doxorubicinol CL4 (L/h)")                                                # Table 2 row CL4, RSE 7%

    # ---------------------------------------------------------------
    # Covariate effects -- BSA was the only retained covariate
    # (Mc Laughlin 2024 Results / Covariate model). Power model:
    #   V2_ind = V2 * (BSA / 1.75)^4.47
    #   V3_ind = V3 * (BSA / 1.75)^11.5
    # The exponents are large because the small central volume of
    # free doxorubicin (V2 ~ 0.5 L) amplifies the BSA dependence
    # (paper text following Table 2). BSA effect on V1 (entrapped)
    # was tested but not retained.
    # ---------------------------------------------------------------
    e_bsa_vc <- 4.47                   ; label("Power exponent of BSA on V2 (free doxorubicin central volume; reference BSA 1.75 m^2)")      # Table 2 row V2_BSA, RSE 19%
    e_bsa_vp <- 11.5                   ; label("Power exponent of BSA on V3 (free doxorubicin peripheral volume; reference BSA 1.75 m^2)")   # Table 2 row V3_BSA, RSE 18%

    # ---------------------------------------------------------------
    # Inter-individual variability. CV% from Mc Laughlin 2024 Table 2
    # converted to log-scale variance via omega^2 = log(1 + CV^2).
    #
    # The paper uses a shared-eta approach between V1 (entrapped) and
    # CL1 (release): the V1 eta is the CL1 eta scaled by the estimated
    # factor theta_shared = 0.643 (Table 2 footnote a, 'theta_shared:
    # Shared eta scale factor for V1'). In this file V1 -> vd and CL1
    # -> krel * vd, so:
    #   log(vd_ind)   = log(vd)   + theta_shared * eta_CL1
    #   log(krel_ind) = log(krel) + (1 - theta_shared) * eta_CL1
    # which gives perfectly-correlated etas on vd and krel with
    # variances scaled by theta_shared^2 and (1 - theta_shared)^2,
    # and covariance = theta_shared * (1 - theta_shared) * omega^2_CL1.
    # Underlying ome^2_CL1 = log(1 + 0.451^2) = 0.1851 (45.1% CV).
    #   var(eta_vd)         = 0.643^2 * 0.1851 = 0.07654
    #   var(eta_krel)       = (1 - 0.643)^2 * 0.1851 = 0.02360
    #   cov(eta_vd,eta_krel)= 0.643 * 0.357 * 0.1851 = 0.04249
    # (implied rho ~ 1.0; this is intentional -- a single underlying
    # random effect drives both vd and krel.)
    # ---------------------------------------------------------------
    etalvd + etalkrel ~ c(0.07654, 0.04249, 0.02360)   # Table 2 rows IIV V1 (28.2% CV), IIV CL1 (45.1% CV), theta_shared = 0.643 (footnote a)
    etalcl       ~ 0.1100                              # Table 2 row IIV CL2 = 34.2% CV -> log(1 + 0.342^2)
    etalcl_doxol ~ 0.02263                             # Table 2 row IIV CL4 = 15.1% CV -> log(1 + 0.151^2)
    # NOTE: Inter-occasion variability (IOV) on CL1 (14.4% CV), V1
    # (8.85% CV), CL2 (22.4% CV), and V2 (126% CV) is reported in
    # Table 2 but NOT encoded structurally here -- nlmixr2lib has no
    # canonical occasion-column convention. See vignette Assumptions
    # and deviations and the Hempel 2003 / Hong 2006 precedents.

    # ---------------------------------------------------------------
    # Residual error. Mc Laughlin 2024 Methods / Stochastic submodel:
    # 'a log-transformed both sides approach with an additive
    # component in the log-domain with separate and uncorrelated
    # parameters for each model species ... best characterised the
    # residual unexplained variability.' Log-additive in NONMEM
    # corresponds to proportional in nlmixr2's linear space (the CV
    # the paper reports back-transforms to the propSd magnitude
    # below).
    # ---------------------------------------------------------------
    propSd_Cc_entrapped <- 0.196       ; label("Proportional residual SD on entrapped doxorubicin concentration (fraction)")                  # Table 2 row RUV Entrapped doxorubicin = 19.6% CV, RSE 5%
    propSd              <- 0.642       ; label("Proportional residual SD on free doxorubicin concentration (fraction)")                       # Table 2 row RUV Free doxorubicin = 64.2% CV, RSE 5%
    propSd_doxol        <- 0.650       ; label("Proportional residual SD on doxorubicinol concentration (fraction)")                          # Table 2 row RUV Doxorubicinol = 65.0% CV, RSE 7%
  })

  model({
    # ---------------------------------------------------------------
    # Individual structural parameters. The shared-eta block in ini()
    # delivers correlated etalvd / etalkrel pairs that reproduce the
    # paper's shared random effect with scale 0.643 on V1 and
    # (1 - 0.643) on log(krel).
    # ---------------------------------------------------------------
    vd        <- exp(lvd        + etalvd)
    krel      <- exp(lkrel      + etalkrel)
    cl        <- exp(lcl        + etalcl)
    cl_doxol  <- exp(lcl_doxol  + etalcl_doxol)
    q         <- exp(lq)
    vc_doxol  <- exp(lvc_doxol)

    # BSA power scaling on the free-doxorubicin volumes (Mc Laughlin
    # 2024 Table 2 footnote b: V2 * (BSA / 1.75)^V2_BSA and
    # V3 * (BSA / 1.75)^V3_BSA).
    vc        <- exp(lvc) * (BSA / 1.75)^e_bsa_vc
    vp        <- exp(lvp) * (BSA / 1.75)^e_bsa_vp

    # ---------------------------------------------------------------
    # ODE system (amounts; mg). Liposomal reservoir 'entrapped'
    # receives the dose; first-order release into the free-doxorubicin
    # central compartment at rate krel. Free doxorubicin is two-
    # compartment with metabolism to doxorubicinol (rate CL/Vc per
    # mass in central). Doxorubicinol is one-compartment with linear
    # elimination.
    #
    # Mc Laughlin 2024 Fig. 2 schematic.
    # ---------------------------------------------------------------
    d/dt(entrapped)     <- -krel * entrapped
    d/dt(central)       <-  krel * entrapped -
                            (cl + q) / vc * central +
                            q / vp * peripheral1
    d/dt(peripheral1)   <-  q / vc * central - q / vp * peripheral1
    d/dt(central_doxol) <-  cl / vc * central -
                            cl_doxol / vc_doxol * central_doxol

    # ---------------------------------------------------------------
    # Observation variables. Doses in mg, volumes in L -> mg/L (the
    # paper's reporting scale). Cc is the canonical free-doxorubicin
    # central observation; Cc_entrapped and Cc_doxol are the per-
    # species concentrations reported in Mc Laughlin 2024 Fig. 4
    # (panels a, b, c respectively).
    # ---------------------------------------------------------------
    Cc_entrapped <- entrapped     / vd
    Cc           <- central       / vc
    Cc_doxol     <- central_doxol / vc_doxol

    # Residual error -- log-additive on each species (see ini() note).
    Cc_entrapped ~ prop(propSd_Cc_entrapped)
    Cc           ~ prop(propSd)
    Cc_doxol     ~ prop(propSd_doxol)
  })
}
