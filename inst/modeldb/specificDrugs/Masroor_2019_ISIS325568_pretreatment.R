Masroor_2019_ISIS325568_pretreatment <- function() {
  description <- paste(
    "QSP. Glucose-insulin-glucagon model of the glucagon challenge test with",
    "a hepatic glucagon-receptor binding / internalization / recycling",
    "sub-model (Masroor 2019), BEFORE treatment with the glucagon-receptor",
    "antisense oligonucleotide ISIS 325568: typical values are the averages",
    "of the eight healthy subjects' individual pre-treatment fits (Table 4).",
    "Covers the 3-h challenge period, in which somatostatin suppresses",
    "endogenous insulin and glucagon secretion while glucose, insulin and",
    "glucagon are infused at constant rates; time 0 is the start of the",
    "challenge (hour 3 of the test) and all states start at the pre-challenge",
    "steady state. Hepatic glucose production is a Hill function (n = 2) of",
    "the glucagon-bound receptor fraction; glucose clearance has",
    "insulin-independent and insulin-dependent Michaelis-Menten terms.",
    "Deterministic: the source fitted each subject separately (smooth",
    "profiling), so no IIV or residual error is reported. The",
    "post-treatment parameter set is Masroor_2019_ISIS325568_posttreatment.",
    sep = " "
  )
  reference <- paste(
    "Masroor S, van Dongen MGJ, Alvarez-Jimenez R, Burggraaf K,",
    "Peletier LA, Peletier MA (2019). Mathematical modeling of the glucagon",
    "challenge test. J Pharmacokinet Pharmacodyn 46(6):553-.",
    "doi:10.1007/s10928-019-09655-2. PMCID: PMC6868112.",
    "Per-subject estimates are in the Electronic Supplementary Material",
    "(Table 7); clinical study: van Dongen MG et al. (2015) J Clin Pharmacol",
    "55(3):298-306, doi:10.1002/jcph.396.",
    sep = " "
  )
  vignette <- "Masroor_2019_glucagon_challenge"

  # The receptor states are fractions of the (unidentifiable) total hepatic
  # glucagon-receptor pool Rtot; the internalized fraction ri = 1 - r - re is
  # eliminated by the conservation law (Masroor 2019 Eq 5).
  paper_specific_compartments <- c(
    "glucagon",
    "target_gcgr_liver",
    "complex_gcgr_liver"
  )

  # The three challenge infusions (Table 3) are dosed into these states.
  dosing <- c("glucose", "insulin", "glucagon")

  units <- list(
    time = "h",
    dosing = paste(
      "Glucose in g (cmt = glucose), insulin in mU (cmt = insulin), glucagon",
      "in pmol (cmt = glucagon). The published challenge infuses 0.24 g/h",
      "glucose, 480 mU/h insulin and 4134.7 pmol/h glucagon for 3 h from",
      "time 0 (Table 3).",
      sep = " "
    ),
    concentration = paste(
      "Glucose Gc in g/L; insulin Ic in mU/L; glucagon GNc in pmol/L;",
      "receptor outputs r_free, r_bound and r_internal (the paper's r, re",
      "and ri) are fractions of the total receptor pool (unitless); hepatic glucose production fhgp in g/L/h.",
      sep = " "
    )
  )

  compartmentData <- list(
    glucose = list(analyte = "glucose", units = "g", specimen = "plasma", verified = TRUE),
    insulin = list(analyte = "insulin", units = "mU", specimen = "plasma", verified = TRUE),
    glucagon = list(analyte = "glucagon", units = "pmol", specimen = "plasma", verified = TRUE),
    target_gcgr_liver = list(
      analyte = "free glucagon receptor (fraction of total)",
      units = "fraction",
      specimen = "tissue",
      verified = TRUE
    ),
    complex_gcgr_liver = list(
      analyte = "glucagon-bound glucagon receptor (fraction of total)",
      units = "fraction",
      specimen = "tissue",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 8,
    n_studies = 1,
    disease_state = "healthy volunteers",
    dose_range = paste(
      "Glucagon challenge test (6 h: 3-h baseline, then 3 h of somatostatin",
      "with insulin, glucagon and glucose infusions), performed before and",
      "after 6 weeks of ISIS 325568 400 mg/week (eight doses). This file is",
      "the before-treatment fit.",
      sep = " "
    ),
    regions = "Netherlands (Centre for Human Drug Research, Leiden)",
    notes = paste(
      "Eight healthy subjects of the ISIS 325568 400 mg/week cohort of",
      "van Dongen 2015 (subjects 30, 31, 32, 33, 35, 36, 38, 40; Masroor",
      "2019 supplementary Table 7). Each subject's pre- and post-treatment",
      "challenges were fitted separately by the smooth profiling method of",
      "Ramsay et al.; this model carries the averages of the eight",
      "pre-treatment fits. Demographics are not reported in Masroor 2019.",
      sep = " "
    )
  )

  ini({
    # Fitted parameters: averages of the eight individual before-treatment fits
    lhgp_base <- log(0.0804); label("Baseline (receptor-independent) hepatic glucose production bG (g/L/h)") # Table 4 'bG' Before average = 0.0804 (SD 0.228)
    lvmax_hgp <- log(5.63); label("Maximal glucagon-receptor-dependent hepatic glucose production V1 (g/L/h)") # Table 4 'V1' Before average = 5.63 (SD 1.97)
    lre50_hgp <- log(0.00501); label("Bound-receptor fraction at half-maximal glucagon-stimulated HGP K1' = K1/Rtot (unitless)") # Table 4 'K1 prime' Before average = 0.00501 (SD 0.00257)
    hill <- fixed(2); label("Hill coefficient of HGP on the bound-receptor fraction (unitless)") # Methods 'The hepatic glucose production rate': n fixed to 2
    lvmax_ii <- log(0.896); label("Maximal insulin-independent glucose consumption Vii (g/L/h)") # Table 4 'Vii' Before average = 0.896 (SD 0.319)
    km_ii <- fixed(0); label("Michaelis constant of insulin-independent glucose consumption Kii (g/L)") # Table 2 'Kii' = 0; Methods: collapsed to zero and set to 0
    lvmax_id <- log(1.52); label("Maximal insulin-dependent glucose consumption Vid (1/h per mU/L)") # Table 4 'Vid' Before average = 1.52 (SD 0.963)
    lkm_id <- log(11); label("Michaelis constant of insulin-dependent glucose consumption Kid (g/L)") # Table 4 'Kid' Before average = 11 (SD 6.62)
    lkdeg_insulin <- log(30.5); label("Insulin degradation rate constant kdegI (1/h)") # Table 4 'kdegI' Before average = 30.5 (SD 8.27)
    lkdeg_glucagon <- log(5.73); label("Glucagon degradation rate constant kdegE (1/h)") # Table 4 'kdegE' Before average = 5.73 (SD 1.84)
    lkint <- log(21.5); label("Internalization rate constant of the glucagon-bound receptor kin (1/h)") # Table 4 'kin' Before average = 21.5 (SD 15.2)
    lvh <- log(4.65); label("Apparent volume of the hepatic interstitial space Vh (L)") # Table 4 'Vh' Before average = 4.65 (SD 2.28)

    # Literature constants held equal for all subjects
    lkon <- fixed(log(0.0036)); label("Glucagon-receptor association rate constant kon (1/pmol/h)") # Table 2 'kon' = 3.6e-3; supplement Table 7 footnote 6e-5 1/pmol/min
    lkoff <- fixed(log(14.4)); label("Glucagon-receptor dissociation rate constant koff (1/h)") # Table 2 'koff' = 14.4; Methods 0.24 1/min
    lkrec <- fixed(log(0.18)); label("Recycling rate constant of the internalized receptor krec (1/h)") # Table 2 'krec' = 0.18; Methods 0.003 1/min
    kint_free <- fixed(0); label("Internalization rate constant of the free receptor kin prime (1/h)") # Table 2 'kin prime' = 0; Methods: collapsed to zero and set to 0

    # Distribution volumes used to convert the infusion rates (Table 3)
    lv_glucose <- fixed(log(4.44)); label("Distribution volume of glucose VG (L)") # Table 3 'VG' = 4.44 L
    lv_insulin <- fixed(log(1.52)); label("Distribution volume of insulin VI (L)") # Table 3 'VI' = 1.52 L
    lv_glucagon <- fixed(log(19.6)); label("Distribution volume of glucagon VE (L)") # Table 3 'VE' = 19.6 L (Table 2 prints 9.6; see vignette)

    # Pre-challenge steady-state concentrations, read from the average-parameter
    # curves of Figure 5 (left column, 0-3 h); not tabulated in the source.
    bl_glucagon <- 17.4; label("Pre-challenge baseline plasma glucagon E-bar (pmol/L)") # digitised from Figure 5 'Before Treatment' E panel, 0-3 h
    bl_insulin <- 10.5; label("Pre-challenge baseline plasma insulin I-bar (mU/L)") # digitised from Figure 5 'Before Treatment' I panel, 0-3 h
  })

  model({
    hgp_base <- exp(lhgp_base)
    vmax_hgp <- exp(lvmax_hgp)
    re50_hgp <- exp(lre50_hgp)
    vmax_ii <- exp(lvmax_ii)
    vmax_id <- exp(lvmax_id)
    km_id <- exp(lkm_id)
    kdeg_insulin <- exp(lkdeg_insulin)
    kdeg_glucagon <- exp(lkdeg_glucagon)
    kint <- exp(lkint)
    vh <- exp(lvh)
    kon <- exp(lkon)
    koff <- exp(lkoff)
    krec <- exp(lkrec)
    v_glucose <- exp(lv_glucose)
    v_insulin <- exp(lv_insulin)
    v_glucagon <- exp(lv_glucagon)

    # Pre-challenge steady state of the receptor sub-model at E = bl_glucagon
    # (Eq 6 / Eq 7d-7e).
    kbind0 <- kon * vh * bl_glucagon
    r0 <- krec * (koff + kint) /
      ((kint_free + krec) * (koff + kint) + (kint + krec) * kbind0)
    re0 <- kbind0 * r0 / (koff + kint)

    # Pre-challenge glucose: positive root of Eq 9 (HGP = total glucose
    # consumption at I = bl_insulin). With km_ii = 0 the root reduces to
    # km_id * (fhgp0 - vmax_ii) / (vmax_ii + vmax_id * bl_insulin - fhgp0).
    fhgp0 <- hgp_base + vmax_hgp * re0^hill / (re50_hgp^hill + re0^hill)
    qa <- fhgp0 - vmax_ii - vmax_id * bl_insulin
    qb <- fhgp0 * (km_ii + km_id) - vmax_ii * km_id - vmax_id * bl_insulin * km_ii
    qc <- fhgp0 * km_ii * km_id
    bl_glucose <- (-qb - sqrt(qb^2 - 4 * qa * qc)) / (2 * qa)

    glucose(0) <- bl_glucose * v_glucose
    insulin(0) <- bl_insulin * v_insulin
    glucagon(0) <- bl_glucagon * v_glucagon
    target_gcgr_liver(0) <- r0
    complex_gcgr_liver(0) <- re0

    Gc <- glucose / v_glucose
    Ic <- insulin / v_insulin
    GNc <- glucagon / v_glucagon

    # Hepatic glucose production (Methods, Hill-type form with n = hill)
    fhgp <- hgp_base + vmax_hgp * complex_gcgr_liver^hill /
      (re50_hgp^hill + complex_gcgr_liver^hill)

    # Eq 8a-8e in amount form for the three plasma species (state = conc * V,
    # so an infusion of Q per hour into the state reproduces the Q/V term).
    # Endogenous insulin and glucagon secretion is zero throughout (somatostatin).
    d/dt(glucose) <- v_glucose * (fhgp - vmax_ii * Gc / (km_ii + Gc) -
      vmax_id * Ic * Gc / (km_id + Gc))
    d/dt(insulin) <- -kdeg_insulin * insulin
    d/dt(glucagon) <- -kdeg_glucagon * glucagon
    d/dt(target_gcgr_liver) <- -kon * vh * GNc * target_gcgr_liver +
      koff * complex_gcgr_liver - kint_free * target_gcgr_liver +
      krec * (1 - target_gcgr_liver - complex_gcgr_liver)
    d/dt(complex_gcgr_liver) <- kon * vh * GNc * target_gcgr_liver -
      (koff + kint) * complex_gcgr_liver

    r_free <- target_gcgr_liver
    r_bound <- complex_gcgr_liver
    r_internal <- 1 - target_gcgr_liver - complex_gcgr_liver
  })
}
