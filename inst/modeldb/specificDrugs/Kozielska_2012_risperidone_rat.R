Kozielska_2012_risperidone_rat <- function() {
  description <- paste(
    "Preclinical (rat).",
    "Hybrid physiology- and mechanism-based population PK-PD model of",
    "risperidone and its active metabolite paliperidone (9-hydroxyrisperidone)",
    "and their striatal dopamine D2 and frontal-cortex serotonin 5-HT2A",
    "receptor occupancy in rats (Kozielska 2012). Plasma PK is a",
    "two-compartment model for each analyte, with risperidone clearance split",
    "into metabolic clearance to paliperidone and clearance by other routes.",
    "IV and IP doses enter the central compartment directly; a fraction of an",
    "IP risperidone dose reaches the paliperidone central compartment as",
    "first-pass metabolism; SC doses are absorbed by consecutive zero- and",
    "first-order processes with a bioavailability shared by both drugs.",
    "Each drug then distributes from plasma into brain vascular and brain",
    "extravascular compartments (passive BBB clearance on unbound drug plus",
    "linear active efflux) and on into striatum and frontal-cortex",
    "compartments, where the two drugs compete for the same D2 and 5-HT2A",
    "receptors with explicit association / dissociation kinetics and a",
    "receptor density (Bmax). Receptor occupancy is the sum of the occupancy",
    "by both drugs. All structural parameters are per kg body weight."
  )
  reference <- paste(
    "Kozielska M, Johnson M, Pilla Reddy V, Vermeulen A, Li C, Grimwood S,",
    "de Greef R, Groothuis GMM, Danhof M, Proost JH.",
    "Pharmacokinetic-pharmacodynamic modeling of the D2 and 5-HT2A receptor",
    "occupancy of risperidone and paliperidone in rats.",
    "Pharm Res. 2012;29(7):1932-1948. doi:10.1007/s11095-012-0722-8."
  )
  vignette <- "Kozielska_2012_risperidone_rat"
  # SC doses go to the depots; IV and IP doses to the central compartments
  # (an IP risperidone dose to both central and central_9oh).
  dosing <- c("depot", "central", "depot_9oh", "central_9oh")
  units <- list(
    time = "h",
    dosing = "mg/kg",
    concentration = "mg/L"
  )

  covariateData <- list(
    ROUTE_IP = list(
      description = paste(
        "Intraperitoneal-route indicator: 1 = the dose record is",
        "intraperitoneal (IP), 0 = intravenous (IV) or subcutaneous (SC)."
      ),
      units = "(binary)",
      type = "binary",
      reference_category = "0 (IV, or SC into a depot)",
      notes = paste(
        "Per-dose-record covariate. An IP risperidone dose is entered as TWO",
        "dose records of the same amount, one into `central` and one into",
        "`central_9oh`, both with ROUTE_IP = 1: the central record takes",
        "bioavailability FIP * (1 - FrFPM) and the central_9oh record",
        "FIP * FrFPM (first-pass conversion to paliperidone). With",
        "ROUTE_IP = 0 a dose into `central` or `central_9oh` is an IV dose",
        "with F = 1. SC doses go to `depot` (risperidone) or `depot_9oh`",
        "(paliperidone) with `rate = -2`; their bioavailability FSC does not",
        "read ROUTE_IP."
      ),
      source_name = "(implicit; per-study route of administration in Table I)"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "risperidone", units = "mg/kg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "risperidone", units = "mg/kg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "risperidone", units = "mg/kg", specimen = "tissue", verified = TRUE),
    depot_9oh = list(analyte = "paliperidone", units = "mg/kg", specimen = "administration site", verified = TRUE),
    central_9oh = list(analyte = "paliperidone", units = "mg/kg", specimen = "plasma", verified = TRUE),
    peripheral1_9oh = list(analyte = "paliperidone", units = "mg/kg", specimen = "tissue", verified = TRUE),
    brain_vascular = list(analyte = "risperidone", units = "mg/kg", specimen = "whole blood", verified = TRUE),
    brain_extravascular = list(analyte = "risperidone", units = "mg/kg", specimen = "tissue", verified = TRUE),
    brain_striatum = list(analyte = "risperidone", units = "mg/kg", specimen = "tissue", verified = TRUE),
    complex_d2 = list(
      analyte = "risperidone-D2 receptor complex",
      units = "mg/kg",
      specimen = "tissue",
      verified = TRUE
    ),
    brain_cortex = list(analyte = "risperidone", units = "mg/kg", specimen = "tissue", verified = TRUE),
    complex_5ht2a = list(
      analyte = "risperidone-5-HT2A receptor complex",
      units = "mg/kg",
      specimen = "tissue",
      verified = TRUE
    ),
    brain_vascular_9oh = list(analyte = "paliperidone", units = "mg/kg", specimen = "whole blood", verified = TRUE),
    brain_extravascular_9oh = list(analyte = "paliperidone", units = "mg/kg", specimen = "tissue", verified = TRUE),
    brain_striatum_9oh = list(analyte = "paliperidone", units = "mg/kg", specimen = "tissue", verified = TRUE),
    complex_d2_9oh = list(
      analyte = "paliperidone-D2 receptor complex",
      units = "mg/kg",
      specimen = "tissue",
      verified = TRUE
    ),
    brain_cortex_9oh = list(analyte = "paliperidone", units = "mg/kg", specimen = "tissue", verified = TRUE),
    complex_5ht2a_9oh = list(
      analyte = "paliperidone-5-HT2A receptor complex",
      units = "mg/kg",
      specimen = "tissue",
      verified = TRUE
    )
  )

  # Receptor-bound pools named by their binding partner, the
  # `complex_<partner>` shape of `complex_fcrn`: `complex_d2` is drug bound to
  # D2 receptors in the striatum, `complex_5ht2a` drug bound to 5-HT2A
  # receptors in the frontal cortex.
  paper_specific_compartments <- c(
    "complex_d2",
    "complex_5ht2a",
    "complex_d2_9oh",
    "complex_5ht2a_9oh",
    "brain_vascular_9oh",
    "brain_extravascular_9oh",
    "brain_striatum_9oh",
    "brain_cortex_9oh"
  )

  population <- list(
    species = "rat (Sprague-Dawley or Wistar, male)",
    n_subjects = 282L,
    n_studies = 12L,
    age_range = "Adult (age not reported)",
    weight_range = "Not reported (no body weights available; parameters are per kg)",
    sex_female_pct = 0,
    disease_state = "Healthy",
    dose_range = paste(
      "Single doses of risperidone 0.01-40 mg/kg IV, IP or SC, and of",
      "paliperidone 0.16-10 mg/kg SC (Table I)."
    ),
    regions = "Pooled industrial datasets from three sponsors",
    notes = paste(
      "Twelve studies (Table I; 282 rats summed over the study rows) from",
      "Janssen Research and Development (Belgium), Merck Sharp and Dohme",
      "(The Netherlands) and Pfizer Global Research and Development (USA),",
      "contributed to the Dutch Top Institute Pharma mechanism-based PK-PD",
      "platform. Mostly one animal per time point (terminal brain sampling).",
      "Observations: plasma risperidone and paliperidone; total brain,",
      "brain-minus-striatum, brain-minus-frontal-cortex, striatum and cortex",
      "concentrations; D2 receptor occupancy in striatum ([3H]raclopride)",
      "and 5-HT2A receptor occupancy in frontal cortex ([3H]M100907), each by",
      "in vivo or ex vivo binding."
    )
  )

  ini({
    # ----- Plasma PK: Table II (original-dataset estimates) -----
    lfip <- log(0.412)
    label("Bioavailability of an IP risperidone dose FIP (fraction)") # Table II: FIP-RIS = 0.412 (RSE 12%)
    lfdepot <- log(0.816)
    label("Bioavailability of an SC dose FSC, risperidone and paliperidone (fraction)") # Table II: FSC = 0.816 (RSE 9%); Methods: same for both drugs
    lfpm_9oh <- log(0.268)
    label("Fraction of the absorbed IP risperidone dose reaching paliperidone central by first-pass metabolism FrFPM (fraction)") # Table II: FrFPM = 0.268 (RSE 12%)
    lka <- log(2.84)
    label("SC first-order absorption rate constant, risperidone KaSC-RIS (1/h)") # Table II: KaSC-RIS = 2.84 (RSE 20%)
    lka_9oh <- log(1.31)
    label("SC first-order absorption rate constant, paliperidone KaSC-PALI (1/h)") # Table II: KaSC-PALI = 1.31 (RSE 22%)
    ld1 <- log(0.161)
    label("SC zero-order input duration, risperidone DRSC-RIS (h)") # Table II: DRSC-RIS = 0.161 (RSE 55%)
    ld1_9oh <- log(0.162)
    label("SC zero-order input duration, paliperidone DRSC-PALI (h)") # Table II: DRSC-PALI = 0.162 (RSE 54%)
    lvc <- log(1.29)
    label("Central volume, risperidone Vc-RIS (L/kg)") # Table II: Vc-RIS = 1.29 (RSE 6%)
    lcl_nonmet <- log(1.62)
    label("Risperidone clearance by routes other than metabolism to paliperidone CLRIS (L/h/kg)") # Table II: CLRIS = 1.62 (RSE 9%)
    lcl_met <- log(0.775)
    label("Metabolic clearance of risperidone to paliperidone CLmet (L/h/kg)") # Table II: CLmet = 0.775 (RSE 11%)
    lvp <- log(0.169)
    label("Peripheral volume, risperidone Vp-RIS (L/kg)") # Table II: Vp-RIS = 0.169 (RSE 16%)
    lq <- log(0.0882)
    label("Intercompartmental clearance, risperidone QRIS (L/h/kg)") # Table II: QRIS = 0.0882 (RSE 25%)
    lvc_9oh <- log(1.27)
    label("Central volume, paliperidone Vc-PALI (L/kg)") # Table II: Vc-PALI = 1.27 (RSE 15%)
    lcl_9oh <- log(1.06)
    label("Clearance, paliperidone CLPALI (L/h/kg)") # Table II: CLPALI = 1.06 (RSE 10%)
    lvp_9oh <- log(0.251)
    label("Peripheral volume, paliperidone Vp-PALI (L/kg)") # Table II: Vp-PALI = 0.251 (RSE 54%)
    lq_9oh <- log(0.269)
    label("Intercompartmental clearance, paliperidone QPALI (L/h/kg)") # Table II: QPALI = 0.269 (RSE 128%)

    # ----- Brain PK and receptor binding: Table III -----
    lcl_bev <- log(2.13)
    label("Passive BBB clearance CLbev, shared by risperidone and paliperidone (L/h/kg)") # Table III: CLbev-RIS = CLbev-PALI = 2.13 (RSE 29%)
    lcl_efflux <- log(9.97)
    label("Active BBB efflux clearance, risperidone CLefflux-RIS (L/h/kg)") # Table III: CLefflux-RIS = 9.97 (RSE 28%)
    lcl_efflux_9oh <- log(47.0)
    label("Active BBB efflux clearance, paliperidone CLefflux-PALI (L/h/kg)") # Table III: CLefflux-PALI = 47.0 (RSE 28%)
    lkd_d2 <- log(0.463)
    label("D2 receptor equilibrium dissociation constant KdD2, both drugs (nM)") # Table III: KdD2 = 0.463 (RSE 14%)
    lkoff_d2 <- log(0.671)
    label("D2 receptor dissociation rate constant koffD2, both drugs (1/h)") # Table III: koffD2 = 0.671 (RSE 19%); konD2 = koff/Kd = 1.45 1/nM/h
    lbmax_d2 <- log(245)
    label("D2 receptor density in striatum Bmax-D2 (nM)") # Table III: Bmax-D2 = 245 (RSE 15%)
    lkd_5ht2a <- log(0.219)
    label("5-HT2A receptor equilibrium dissociation constant Kd5-HT2A, both drugs (nM)") # Table III: Kd5-HT2A = 0.219 (RSE 15%)
    lkoff_5ht2a <- log(0.525)
    label("5-HT2A receptor dissociation rate constant koff5-HT2A, both drugs (1/h)") # Table III: koff5-HT2A = 0.525 (RSE 25%); kon5-HT2A = koff/Kd = 2.40 1/nM/h
    lbmax_5ht2a <- log(46.5)
    label("5-HT2A receptor density in frontal cortex Bmax-5-HT2A (nM)") # Table III: Bmax-5-HT2A = 46.5 (RSE 11%)

    # ----- Fixed physiological values and literature binding (Methods, PK-PD Modeling) -----
    lvbv <- fixed(log(0.00024))
    label("Brain vascular volume Vbv (L/kg)") # Methods: 0.00024 L/kg (refs 17, 18)
    lvbev <- fixed(log(0.00656))
    label("Brain extravascular volume Vbev (L/kg)") # Methods: total extra-vascular 0.00656 L/kg (refs 17, 18)
    lvstr <- fixed(log(0.0002))
    label("Striatum volume Vstr (L/kg)") # Methods: 0.0002 L/kg (refs 17, 18)
    lvcor <- fixed(log(0.0035))
    label("Frontal cortex volume Vcor (L/kg)") # Methods: 0.0035 L/kg (ref 22)
    lcl_bv <- fixed(log(0.312))
    label("Plasma to brain vascular clearance CLbv = cerebral blood flow, both drugs (L/h/kg)") # Methods: 0.312 L/h/kg (ref 17)
    lcl_str <- fixed(log(500))
    label("Brain extravascular to striatum clearance CLstr (L/h/kg)") # Methods: fixed to a high value, 500 L/h/kg
    lcl_cor <- fixed(log(500))
    label("Brain extravascular to frontal cortex clearance CLcor (L/h/kg)") # Methods: 'As for striatum, we assumed fast equilibration'; value assumed equal to CLstr = 500
    fu_plasma <- fixed(0.0798)
    label("Unbound fraction in plasma, risperidone (fraction)") # Methods: fuplasma-RIS = 0.0798 (ref 20)
    fu_brain <- fixed(0.0699)
    label("Unbound fraction in brain, risperidone (fraction)") # Methods: fubrain-RIS = 0.0699 (ref 20)
    fu_plasma_9oh <- fixed(0.129)
    label("Unbound fraction in plasma, paliperidone (fraction)") # Methods: fuplasma-PALI = 0.129 (ref 20)
    fu_brain_9oh <- fixed(0.0755)
    label("Unbound fraction in brain, paliperidone (fraction)") # Methods: fubrain-PALI = 0.0755 (ref 20)

    # ----- IIV: Table II, reported as %CV; omega^2 = (CV/100)^2 -----
    # The brain PK-PD parameters carry no IIV (Results, PK-PD: fixed to zero).
    etalfip ~ 0.6496 # Table II: IIV-FIP 80.6% CV
    etalka ~ 0.2153 # Table II: IIV-KaSC-RIS 46.4% CV
    etalka_9oh ~ 0.1183 # Table II: IIV-KaSC-PALI 34.4% CV
    etald1 ~ 0.8317 # Table II: IIV-DRSC-RIS 91.2% CV
    etalcl_nonmet ~ 0.0930 # Table II: IIV-CLRIS 30.5% CV
    etalcl_9oh ~ 0.0262 # Table II: IIV-CLPALI 16.2% CV

    # ----- Residual error -----
    # Plasma and brain concentrations were fitted log-transformed, so the
    # tabulated proportional errors are SDs on the log scale (lnorm).
    expSd <- 0.233
    label("Plasma risperidone residual SD (log scale)") # Table II: proportional residual error, risperidone 0.233
    expSd_9oh <- 0.186
    label("Plasma paliperidone residual SD (log scale)") # Table II: proportional residual error, paliperidone 0.186
    expSd_Cbrain <- 0.362
    label("Brain risperidone residual SD (log scale)") # Table III: proportional BC RIS 0.362
    expSd_Cbrain_9oh <- 0.424
    label("Brain paliperidone residual SD (log scale)") # Table III: proportional BC PALI 0.424
    addSd_D2RO <- 0.177
    label("D2 receptor occupancy additive residual SD (fraction)") # Table III: additive % D2 RO 17.7 (percentage points)
    addSd_RO_5HT2A <- 0.182
    label("5-HT2A receptor occupancy additive residual SD (fraction)") # Table III: additive % 5-HT2A RO 18.2 (percentage points)
  })

  model({
    # Molecular weights (g/mol) convert bound concentrations from mg/L to nM;
    # not given in the paper (formulas C23H27FN4O2 and C23H27FN4O3).
    mw <- 410.49
    mw_9oh <- 426.49

    # ----- Individual parameters -----
    fip <- exp(lfip + etalfip)
    fsc <- exp(lfdepot)
    fpm_9oh <- exp(lfpm_9oh)
    ka <- exp(lka + etalka)
    ka_9oh <- exp(lka_9oh + etalka_9oh)
    d1 <- exp(ld1 + etald1)
    d1_9oh <- exp(ld1_9oh)
    vc <- exp(lvc)
    cl_nonmet <- exp(lcl_nonmet + etalcl_nonmet)
    cl_met <- exp(lcl_met)
    vp <- exp(lvp)
    q <- exp(lq)
    vc_9oh <- exp(lvc_9oh)
    cl_9oh <- exp(lcl_9oh + etalcl_9oh)
    vp_9oh <- exp(lvp_9oh)
    q_9oh <- exp(lq_9oh)

    vbv <- exp(lvbv)
    vbev <- exp(lvbev)
    vstr <- exp(lvstr)
    vcor <- exp(lvcor)
    cl_bv <- exp(lcl_bv)
    cl_bev <- exp(lcl_bev)
    cl_efflux <- exp(lcl_efflux)
    cl_efflux_9oh <- exp(lcl_efflux_9oh)
    cl_str <- exp(lcl_str)
    cl_cor <- exp(lcl_cor)
    kd_d2 <- exp(lkd_d2)
    koff_d2 <- exp(lkoff_d2)
    kon_d2 <- koff_d2 / kd_d2
    bmax_d2 <- exp(lbmax_d2)
    kd_5ht2a <- exp(lkd_5ht2a)
    koff_5ht2a <- exp(lkoff_5ht2a)
    kon_5ht2a <- koff_5ht2a / kd_5ht2a
    bmax_5ht2a <- exp(lbmax_5ht2a)

    # ----- Concentrations (mg/L) -----
    Cp <- central / vc
    Cp_9oh <- central_9oh / vc_9oh
    Cbv <- brain_vascular / vbv
    Cbv_9oh <- brain_vascular_9oh / vbv
    Cbev <- brain_extravascular / vbev
    Cbev_9oh <- brain_extravascular_9oh / vbev
    Cstr <- brain_striatum / vstr
    Cstr_9oh <- brain_striatum_9oh / vstr
    Ccor <- brain_cortex / vcor
    Ccor_9oh <- brain_cortex_9oh / vcor

    # Receptor-bound concentrations (nM): Appendix CB = 1000 * A / V / MW,
    # written there for amounts in ug; 1e6 / MW for amounts in mg.
    cb_d2 <- complex_d2 / vstr * 1e6 / mw
    cb_d2_9oh <- complex_d2_9oh / vstr * 1e6 / mw_9oh
    cb_5ht2a <- complex_5ht2a / vcor * 1e6 / mw
    cb_5ht2a_9oh <- complex_5ht2a_9oh / vcor * 1e6 / mw_9oh
    free_d2 <- bmax_d2 - cb_d2 - cb_d2_9oh
    free_5ht2a <- bmax_5ht2a - cb_5ht2a - cb_5ht2a_9oh

    # ----- Plasma PK (Appendix) -----
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - q * (Cp - peripheral1 / vp) -
      cl_nonmet * Cp - cl_met * Cp
    d/dt(peripheral1) <- q * (Cp - peripheral1 / vp)
    d/dt(depot_9oh) <- -ka_9oh * depot_9oh
    d/dt(central_9oh) <- ka_9oh * depot_9oh -
      q_9oh * (Cp_9oh - peripheral1_9oh / vp_9oh) -
      cl_9oh * Cp_9oh + cl_met * Cp
    d/dt(peripheral1_9oh) <- q_9oh * (Cp_9oh - peripheral1_9oh / vp_9oh)

    # ----- Brain kinetics and receptor binding, risperidone (Appendix) -----
    # Plasma loses no drug to the brain (the plasma equations carry no brain
    # term). Only unbound drug crosses the BBB; efflux is linear.
    d/dt(brain_vascular) <- cl_bv * Cp - cl_bv * Cbv -
      cl_bev * fu_plasma * Cbv + cl_bev * fu_brain * Cbev +
      cl_efflux * fu_brain * Cbev
    d/dt(brain_extravascular) <- cl_bev * fu_plasma * Cbv -
      cl_bev * fu_brain * Cbev - cl_efflux * fu_brain * Cbev -
      cl_str * fu_brain * Cbev + cl_str * fu_brain * Cstr -
      cl_cor * fu_brain * Cbev + cl_cor * fu_brain * Ccor
    d/dt(brain_striatum) <- cl_str * fu_brain * Cbev -
      cl_str * fu_brain * Cstr -
      kon_d2 * fu_brain * brain_striatum * free_d2 + koff_d2 * complex_d2
    d/dt(complex_d2) <- kon_d2 * fu_brain * brain_striatum * free_d2 -
      koff_d2 * complex_d2
    d/dt(brain_cortex) <- cl_cor * fu_brain * Cbev -
      cl_cor * fu_brain * Ccor -
      kon_5ht2a * fu_brain * brain_cortex * free_5ht2a +
      koff_5ht2a * complex_5ht2a
    d/dt(complex_5ht2a) <- kon_5ht2a * fu_brain * brain_cortex * free_5ht2a -
      koff_5ht2a * complex_5ht2a

    # ----- Brain kinetics and receptor binding, paliperidone (Appendix) -----
    d/dt(brain_vascular_9oh) <- cl_bv * Cp_9oh - cl_bv * Cbv_9oh -
      cl_bev * fu_plasma_9oh * Cbv_9oh + cl_bev * fu_brain_9oh * Cbev_9oh +
      cl_efflux_9oh * fu_brain_9oh * Cbev_9oh
    d/dt(brain_extravascular_9oh) <- cl_bev * fu_plasma_9oh * Cbv_9oh -
      cl_bev * fu_brain_9oh * Cbev_9oh -
      cl_efflux_9oh * fu_brain_9oh * Cbev_9oh -
      cl_str * fu_brain_9oh * Cbev_9oh + cl_str * fu_brain_9oh * Cstr_9oh -
      cl_cor * fu_brain_9oh * Cbev_9oh + cl_cor * fu_brain_9oh * Ccor_9oh
    d/dt(brain_striatum_9oh) <- cl_str * fu_brain_9oh * Cbev_9oh -
      cl_str * fu_brain_9oh * Cstr_9oh -
      kon_d2 * fu_brain_9oh * brain_striatum_9oh * free_d2 +
      koff_d2 * complex_d2_9oh
    d/dt(complex_d2_9oh) <- kon_d2 * fu_brain_9oh * brain_striatum_9oh * free_d2 -
      koff_d2 * complex_d2_9oh
    d/dt(brain_cortex_9oh) <- cl_cor * fu_brain_9oh * Cbev_9oh -
      cl_cor * fu_brain_9oh * Ccor_9oh -
      kon_5ht2a * fu_brain_9oh * brain_cortex_9oh * free_5ht2a +
      koff_5ht2a * complex_5ht2a_9oh
    d/dt(complex_5ht2a_9oh) <- kon_5ht2a * fu_brain_9oh * brain_cortex_9oh * free_5ht2a -
      koff_5ht2a * complex_5ht2a_9oh

    # ----- Dosing -----
    # SC: consecutive zero-order (duration d1) and first-order (ka) input.
    f(depot) <- fsc
    dur(depot) <- d1
    f(depot_9oh) <- fsc
    dur(depot_9oh) <- d1_9oh
    # IP risperidone: FIP * (1 - FrFPM) reaches risperidone central and
    # FIP * FrFPM paliperidone central (the same dose entered on both).
    # ROUTE_IP = 0 gives an IV dose with F = 1.
    f(central) <- 1 + ROUTE_IP * (fip * (1 - fpm_9oh) - 1)
    f(central_9oh) <- 1 + ROUTE_IP * (fip * fpm_9oh - 1)

    # ----- Observations -----
    # Brain homogenate concentration: total amount over total volume of the
    # brain compartments sampled.
    a_brain <- brain_vascular + brain_extravascular + brain_striatum +
      complex_d2 + brain_cortex + complex_5ht2a
    a_brain_9oh <- brain_vascular_9oh + brain_extravascular_9oh +
      brain_striatum_9oh + complex_d2_9oh + brain_cortex_9oh +
      complex_5ht2a_9oh
    v_brain <- vbv + vbev + vstr + vcor
    Cbrain_exstr <- (a_brain - brain_striatum - complex_d2) / (v_brain - vstr)
    Cbrain_excor <- (a_brain - brain_cortex - complex_5ht2a) / (v_brain - vcor)
    Cstriatum <- (brain_striatum + complex_d2) / vstr
    Ccortex <- (brain_cortex + complex_5ht2a) / vcor
    Cstriatum_9oh <- (brain_striatum_9oh + complex_d2_9oh) / vstr
    Ccortex_9oh <- (brain_cortex_9oh + complex_5ht2a_9oh) / vcor

    Cc <- Cp
    Cc_9oh <- Cp_9oh
    Cbrain <- a_brain / v_brain
    Cbrain_9oh <- a_brain_9oh / v_brain
    D2RO <- (cb_d2 + cb_d2_9oh) / bmax_d2
    RO_5HT2A <- (cb_5ht2a + cb_5ht2a_9oh) / bmax_5ht2a

    Cc ~ lnorm(expSd)
    Cc_9oh ~ lnorm(expSd_9oh)
    Cbrain ~ lnorm(expSd_Cbrain)
    Cbrain_9oh ~ lnorm(expSd_Cbrain_9oh)
    D2RO ~ add(addSd_D2RO)
    RO_5HT2A ~ add(addSd_RO_5HT2A)
  })
}
