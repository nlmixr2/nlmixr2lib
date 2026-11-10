Wang_2024_docetaxel <- function() {
  description <- paste(
    "Three-compartment population PK model for total and unbound docetaxel",
    "in men with metastatic prostate cancer after a 1-h intravenous infusion",
    "of standard docetaxel and after a single dose of oral docetaxel taken",
    "1 h after 15 mg of the P-glycoprotein inhibitor encequidar (oDox + E;",
    "Wang 2024). The disposition parameters are those of UNBOUND docetaxel:",
    "the central state holds the unbound-equivalent amount, so",
    "Cunbound = central / vc, and total plasma docetaxel is reconstructed as",
    "Cc = Cunbound / fu with a constant fraction unbound that differs",
    "between the intravenous formulation (fu 0.67%, polysorbate-80 micelle)",
    "and oDox + E (fu 1.02%), selected by ROUTE_IV. Oral absorption is",
    "first order with a lag time and a constant logit-normal bioavailability",
    "(F 25%). Residual error is combined additive plus proportional,",
    "estimated separately for each of the four route-by-analyte data sets.",
    "Nine patients; no covariates were tested."
  )
  reference <- paste(
    "Wang D, Jackson C, Hung N, Hung T, Kwan R, Chan W-K, Qin A,",
    "Hughes-Medlicott NJ, Glue P, Duffull S. Oral docetaxel plus",
    "encequidar - a pharmacokinetic model and evaluation against IV",
    "docetaxel. J Pharmacokinet Pharmacodyn. 2024;51(4):335-.",
    "doi:10.1007/s10928-024-09913-y. Correction (author affiliations only):",
    "doi:10.1007/s10928-024-09924-9.",
    sep = " "
  )
  vignette <- "Wang_2024_docetaxel"

  # Wang 2024 Supplement 1 estimates a separate combined additive +
  # proportional residual error for each of four data sets: total and
  # unbound docetaxel, each after IV docetaxel and after oDox + E. The
  # route-suffixed names are declared here so checkModelConventions()
  # accepts them; model() selects between each IV/oral pair with ROUTE_IV.
  paper_specific_residual_sds <- c(
    "propSdIv",
    "addSdIv",
    "propSdOral",
    "addSdOral",
    "propSdIv_Cunbound",
    "addSdIv_Cunbound",
    "propSdOral_Cunbound",
    "addSdOral_Cunbound"
  )

  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    ROUTE_IV = list(
      description = paste(
        "1 = intravenous docetaxel period (1-h infusion into central),",
        "0 = oral docetaxel plus encequidar period (dose into depot).",
        "Selects the fraction unbound (fu_iv vs fu_po) used to reconstruct",
        "total docetaxel and the residual-error pair of both outputs."
      ),
      units = "(binary)",
      type = "binary",
      reference_category = "0 (oral docetaxel plus encequidar)",
      notes = paste(
        "Per-period indicator: every patient received both routes in",
        "separate study periods at least 3 weeks apart, so the column varies",
        "within subject and must be set on observation records as well as",
        "dose records. The reference category is oral, as in",
        "Koolen_2010_docetaxel.R."
      ),
      source_name = "route (IV docetaxel vs oDox + E)"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "docetaxel", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "docetaxel", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "docetaxel", units = "mg", specimen = "tissue", verified = TRUE),
    peripheral2 = list(analyte = "docetaxel", units = "mg", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 9L,
    n_studies = 1L,
    age_range = "not reported in the source paper",
    weight_range = "not reported in the source paper",
    sex_female_pct = 0,
    race_ethnicity = "not reported in the source paper",
    disease_state = "Metastatic prostate cancer, receiving standard-of-care IV docetaxel",
    dose_range = paste(
      "IV docetaxel ~75 mg/m^2 as a 1-h infusion; then, at least 3 weeks later,",
      "a single oral docetaxel dose of 75, 150 or 300 mg/m^2 taken 1 h after",
      "15 mg oral encequidar"
    ),
    regions = "not reported in the source paper",
    notes = paste(
      "Phase I open-label two-period, two-treatment study (Wang 2024 Methods,",
      "'Study design'; trial reported by Jackson et al., reference 21).",
      "Intensive sampling: 24 samples after IV and 23 after oDox + E; total",
      "docetaxel by HPLC-MS/MS (LLOQ 2 ng/mL) in all samples, unbound",
      "docetaxel by ultrafiltration HPLC-MS/MS (LLOQ 0.084 ng/mL) in 7 (IV) and",
      "8 (oral) optimal-design-selected samples per patient; 135 unbound",
      "samples, 25% below LLOQ, handled with the M6 method."
    )
  )

  ini({
    # Wang 2024 Supplement 1, 'Final Model (Combined)' column (sub-model 6).
    # The parameters are reported 'in the context of unbound docetaxel'
    # (Results, 'Final model'): CL, Vc, Q and Vp are unbound-concentration
    # parameters, and total clearance is CL * fu_iv = 8570 * 0.0067 = 57.3 L/h.
    lcl <- log(8570); label("Clearance of unbound docetaxel (L/h)") # Supplement 1, CL(/F) = 8570 L/h (RSE 15%)
    lvc <- log(3500); label("Central volume of unbound docetaxel (L)") # Supplement 1, Vc = 3500 L (RSE 24%)
    lq <- log(1250); label("Intercompartmental clearance central-peripheral1 (L/h)") # Supplement 1, Qp1 = 1250 L/h (RSE 25%)
    lvp <- log(9450); label("First peripheral volume (L)") # Supplement 1, Vp1 = 9450 L (RSE 69%)
    lq2 <- log(2810); label("Intercompartmental clearance central-peripheral2 (L/h)") # Supplement 1, Qp2 = 2810 L/h (RSE 14%)
    lvp2 <- log(61900); label("Second peripheral volume (L)") # Supplement 1, Vp2 = 61900 L (RSE 15%)

    # Oral docetaxel plus encequidar absorption.
    lka <- log(0.223); label("First-order absorption rate constant (1/h)") # Supplement 1, ka = 0.223 1/h (RSE 19%)
    ltlag <- log(0.242); label("Absorption lag time (h)") # Supplement 1, ALAG = 0.242 h (RSE 1%)
    logitfdepot <- qlogis(0.25); label("Logit of oral bioavailability with encequidar (F = 0.25)") # Supplement 1, F = 25% (RSE 30%)

    # Constant plasma binding. The paper estimates binding constants
    # B = 1 / fu and tabulates them as fractions unbound.
    lfu_iv <- log(0.0067); label("Fraction unbound after IV docetaxel (fraction)") # Supplement 1, f unbound IV = 0.67% (RSE 18%)
    lfu_po <- log(0.0102); label("Fraction unbound after oral docetaxel plus encequidar (fraction)") # Supplement 1, f unbound PO = 1.02% (RSE 26%)

    # Between-subject variability. Supplement 1 footnote: omega values are
    # shown as CV% 'calculated as the square root of the omega value', so
    # the variance is the square of the tabulated CV (as a fraction). The
    # logitF entry (0.71, no percent sign) is the SD on the logit scale.
    # The Vp2 omega was estimated at 0% (RSE NA) and is not carried.
    etalcl ~ 0.0729 # Supplement 1, BSV CL = 27% (RSE 64%); 0.27^2
    etalogitfdepot ~ 0.5041 # Supplement 1, BSV logitF = 0.71 (RSE 51%); 0.71^2 on the logit scale
    etalka ~ 0.2209 # Supplement 1, BSV ka = 47% (RSE 60%); 0.47^2
    etaltlag ~ 0.0009 # Supplement 1, BSV ALAG = 3% (RSE 110%); 0.03^2
    etalvc ~ 0.3969 # Supplement 1, BSV Vc = 63% (RSE 51%); 0.63^2
    etalvp ~ 2.1904 # Supplement 1, BSV Vp1 = 148% (RSE 48%); 1.48^2
    etalfu_iv ~ 0.0784 # Supplement 1, BSV BIV = 28% (RSE 62%); 0.28^2, on B = 1/fu
    etalfu_po ~ 0.5184 # Supplement 1, BSV BPO = 72% (RSE 48%); 0.72^2, on B = 1/fu

    # Residual error, combined additive + proportional per data set.
    addSdIv <- 0.001; label("Additive residual SD, total docetaxel after IV (ng/mL)") # Supplement 1, Tot IV add = '<0.001' ng/mL (RSE NA); upper bound used
    propSdIv <- 0.143; label("Proportional residual SD, total docetaxel after IV (fraction)") # Supplement 1, Tot IV prop = 0.143 (RSE 12%)
    addSdOral <- 0.199; label("Additive residual SD, total docetaxel after oDox + E (ng/mL)") # Supplement 1, Tot PO add = 0.199 ng/mL (RSE 181%)
    propSdOral <- 0.322; label("Proportional residual SD, total docetaxel after oDox + E (fraction)") # Supplement 1, Tot PO prop = 0.322 (RSE 14%)
    addSdIv_Cunbound <- 0.005; label("Additive residual SD, unbound docetaxel after IV (ng/mL)") # Supplement 1, UF IV add = 0.005 ng/mL (RSE 139%)
    propSdIv_Cunbound <- 0.538; label("Proportional residual SD, unbound docetaxel after IV (fraction)") # Supplement 1, UF IV prop = 0.538 (RSE 36%)
    addSdOral_Cunbound <- 0.003; label("Additive residual SD, unbound docetaxel after oDox + E (ng/mL)") # Supplement 1, UF PO add = 0.003 ng/mL (RSE 99%)
    propSdOral_Cunbound <- 0.337; label("Proportional residual SD, unbound docetaxel after oDox + E (fraction)") # Supplement 1, UF PO prop = 0.337 (RSE 18%)
  })
  model({
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    q <- exp(lq)
    vp <- exp(lvp + etalvp)
    q2 <- exp(lq2)
    vp2 <- exp(lvp2)
    ka <- exp(lka + etalka)
    tlag <- exp(ltlag + etaltlag)
    fdepot <- expit(logitfdepot + etalogitfdepot)

    # The paper puts exponential BSV on the binding constant B = 1 / fu;
    # exp(eta) on B is exp(-eta) on fu, the same log-normal distribution.
    fu_iv <- exp(lfu_iv + etalfu_iv)
    fu_po <- exp(lfu_po + etalfu_po)
    fu <- ROUTE_IV * fu_iv + (1 - ROUTE_IV) * fu_po

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1 - k13 * central + k31 * peripheral2
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2

    f(depot) <- fdepot
    alag(depot) <- tlag

    # Supplement 2: the central compartment holds plasma unbound docetaxel
    # in equilibrium with the plasma-bound pool through B_IV or B_PO.
    # mg/L -> ng/mL.
    Cunbound <- central / vc * 1000
    Cc <- Cunbound / fu

    propSd <- ROUTE_IV * propSdIv + (1 - ROUTE_IV) * propSdOral
    addSd <- ROUTE_IV * addSdIv + (1 - ROUTE_IV) * addSdOral
    propSd_Cunbound <- ROUTE_IV * propSdIv_Cunbound + (1 - ROUTE_IV) * propSdOral_Cunbound
    addSd_Cunbound <- ROUTE_IV * addSdIv_Cunbound + (1 - ROUTE_IV) * addSdOral_Cunbound

    Cc ~ add(addSd) + prop(propSd)
    Cunbound ~ add(addSd_Cunbound) + prop(propSd_Cunbound)
  })
}
