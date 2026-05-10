Panoilia_2015_bevacizumab <- function() {
  description <- "Quasi-steady-state target-mediated drug-disposition (TMDD QSS) model for IV bevacizumab and free VEGF165 in adults with stage IV colorectal cancer, with fixed allometric body-weight scaling on PK clearances and volumes (Panoilia 2015, Table 3 TMDD model column)"
  reference <- "Panoilia E, Schindler E, Samantas E, Aravantinos G, Kalofonos HP, Christodoulou C, Patrinos GP, Friberg LE, Sivolapenko G. A pharmacokinetic binding model for bevacizumab and VEGF165 in colorectal cancer patients. Cancer Chemother Pharmacol. 2015;75(4):791-803. doi:10.1007/s00280-015-2701-3"
  vignette <- "Panoilia_2015_bevacizumab"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L (bevacizumab) and ng/L (free VEGF165)")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed in source analysis. Allometric power scaling with fixed exponents 0.75 on CL/Q and 1 on V1/V2, centred at 70 kg (study median, Panoilia 2015 Methods Eq 3 and Results page 797).",
      source_name        = "Actual body weight"
    )
  )

  population <- list(
    n_subjects        = 19L,
    n_studies         = 1L,
    age_range         = "37-73 years",
    age_median        = "60 years",
    weight_range      = "50-94 kg",
    weight_median     = "70 kg",
    sex_female_pct    = 42,
    race_ethnicity    = "Greek (single-country observational study); race not formally reported.",
    disease_state     = "Adults with stage IV colorectal cancer (mCRC); ECOG performance status <= 2.",
    dose_range        = "5 mg/kg IV every 2 weeks (BEV-FOLFIRI / BEV-FOLFOX), 7.5 mg/kg IV every 3 weeks (BEV-CAPIRI), or 10 mg/kg every 2 weeks (1 patient).",
    regions           = "Greece",
    n_observations    = "86 total bevacizumab + 93 free VEGF165 serum concentrations (median 4 samples per patient, range 2-10).",
    co_medication     = "BEV-FOLFIRI (5-FU/leucovorin/irinotecan), BEV-FOLFOX (5-FU/leucovorin/oxaliplatin), or BEV-CAPIRI (capecitabine/irinotecan).",
    notes             = "Three-centre prospective observational study in Greece. SNPs (rs699947 / -2578C/A, rs1570360 / -1154G/A, rs2010963 / -634G/C) were tested but no SNP covariate was retained in the final TMDD model after randomization-test correction (Panoilia 2015 Discussion)."
  )

  ini({
    # Structural PK parameters (reference 70 kg adult)
    lcl   <- log(0.18);   label("Linear clearance of free bevacizumab for a 70 kg adult (CL, L/day)")  # Table 3 TMDD column row CL
    lvc   <- log(3.23);   label("Central volume of distribution (V1, L)")                              # Table 3 TMDD column row V1
    lq    <- log(1.38);   label("Inter-compartmental clearance (Q, L/day)")                            # Table 3 TMDD column row Q
    lvp   <- log(3.1);    label("Peripheral volume of distribution (V2, L)")                           # Table 3 TMDD column row V2

    # Target turnover (free VEGF165). BM0 reported as 0.0053 nM (= 212 ng/L) under the paper's
    # 1:1 monomeric bevacizumab-VEGF165 binding assumption.
    lkout <- log(0.401);  label("First-order elimination rate constant of free VEGF165 (kout, 1/day)") # Table 3 row kout
    lbm0  <- log(0.0053); label("Baseline (pre-dose) free VEGF165 concentration (BM0, nM)")            # Table 3 row BM0 (= 212 ng/L)

    # QSS dissociation constant (same units as Ctot, Rtot)
    lkss  <- log(267);    label("Quasi-steady-state dissociation constant for bevacizumab-VEGF165 binding (Kss, nM)")  # Table 3 row Kss

    # Fixed allometric exponents on body weight. Panoilia 2015 Methods (Eq 3): the power exponent
    # k was "either estimated or fixed to a certain value (0.75 for clearance and 1 for volume
    # parameters)". Discussion: "the developed TMDD (binding) model did not include any covariates
    # except for body weight in all clearance and volume parameters". The shared exponents follow.
    e_wt_cl_q  <- fixed(0.75); label("Fixed allometric exponent on log(WT/70) shared by CL and Q (unitless)")
    e_wt_vc_vp <- fixed(1);    label("Fixed allometric exponent on log(WT/70) shared by V1 and V2 (unitless)")

    # IIV (exponential, log-normal). omega^2 = log(CV^2 + 1) for CV reported as % in Table 3.
    # CL CV 20%   -> omega^2 = log(1 + 0.20^2) = 0.03922
    # V1 CV 22%   -> omega^2 = log(1 + 0.22^2) = 0.04727
    # BM0 CV 24%  -> omega^2 = log(1 + 0.24^2) = 0.05599
    etalcl  ~ 0.03922
    etalvc  ~ 0.04727
    etalbm0 ~ 0.05599

    # Residual error: paper Methods page 796 ("Modeling was performed on log-transformed data ...
    # an additive error model on log-transformed data"). NONMEM additive on log scale is
    # proportional in linear space; Table 3 reports the resulting CV % directly.
    propSd          <- 0.28; label("Proportional residual error for total bevacizumab (fraction)")  # Table 3 row Prop. error_bev
    propSd_freeVegf <- 0.32; label("Proportional residual error for free VEGF165 (fraction)")        # Table 3 row Prop. error_VEGF165
  })
  model({
    # ---- Unit-conversion constants ----
    # Bevacizumab MW (g/mol). Standard reference value for the marketed antibody; not stated
    # numerically in Panoilia 2015 but required to convert mg dose into the molar concentrations
    # used by the TMDD binding equations. Documented as a model assumption in the vignette.
    mw_bev <- 149000

    # VEGF165 MW (g/mol). Back-calculated from Panoilia 2015 Table 3 BM0 = 0.0053 nM = 212 ng/L
    # (1:1 binding assumption): MW = 212 / 0.0053 = 40000 g/mol, consistent with the homodimeric
    # VEGF165 form (~22 kDa monomer).
    mw_vegf165 <- 40000

    # mg-to-nmol conversion factor for bevacizumab (1 mg = 1e6 / MW nmol).
    nmol_per_mg_bev <- 1e6 / mw_bev

    # ---- Individual parameters (Panoilia 2015 Methods Eq 3 with k fixed) ----
    cl   <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl_q
    vc   <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc_vp
    q    <- exp(lq)           * (WT / 70)^e_wt_cl_q
    vp   <- exp(lvp)          * (WT / 70)^e_wt_vc_vp
    kout <- exp(lkout)
    bm0  <- exp(lbm0 + etalbm0)
    kss  <- exp(lkss)

    # Zero-order target synthesis (kin = kout * BM0 at steady state with no drug; paper Results:
    # "kin, which is defined as the typical value of BM0 times the typical value of kout").
    kin  <- kout * bm0

    # Complex elimination rate constant. Panoilia 2015 Results: "the elimination clearance of the
    # bevacizumab-VEGF165 complex ... was therefore set equal to the CL of the free bevacizumab
    # (CL_RC = CL)." With CL_RC = CL, the per-concentration rate constant is kint = CL_RC / V1 = CL / V1.
    kint <- cl / vc

    # ---- QSS algebra (Gibiansky et al. 2008 / Panoilia 2015 Eq 7) ----
    # central tracks TOTAL bevacizumab amount (free + complex) in nmol;
    # total_target tracks TOTAL (free + bound) VEGF165 concentration in V1 (nM).
    ctot <- central / vc
    ttot <- total_target

    # Numerically stable QSS quadratic. Panoilia 2015 Eq 7:
    #   C = 0.5 * { (Ctot - Rtot - Kss) + sqrt((Ctot - Rtot - Kss)^2 + 4 * Kss * Ctot) }
    # Equivalently (and stabler), the bevacizumab-VEGF165 complex is the negative root of the
    # mass-balance quadratic with discriminant
    #   (Ctot + Ttot + Kss)^2 - 4 * Ctot * Ttot
    # = (Ctot - Ttot)^2 + 2 * (Ctot + Ttot) * Kss + Kss^2 (always >= 0).
    qss_disc <- (ctot - ttot)^2 + 2 * (ctot + ttot) * kss + kss^2
    complex  <- ((ctot + ttot + kss) - sqrt(qss_disc)) / 2
    cfree    <- ctot - complex
    tfree    <- ttot - complex

    # ---- Dosing: mg input -> nmol via the bioavailability multiplier on central ----
    f(central) <- nmol_per_mg_bev

    # ---- ODEs (Panoilia 2015 Eqs 5, 6) ----
    # Central drug: free + complex; complex eliminates at rate kint * V1 = CL.
    # Peripheral compartment carries free drug only.
    d/dt(central)     <- -cl * cfree - q * cfree + q * (peripheral1 / vp) - kint * complex * vc
    d/dt(peripheral1) <-               q * cfree - q * (peripheral1 / vp)

    # Total target: kin - kout * R_free - kint * complex (Eq 6 expanded; QSS-TMDD convention
    # confines the target to the central compartment).
    d/dt(total_target) <- kin - kout * tfree - kint * complex
    total_target(0)    <- bm0  # steady-state baseline before drug (no drug -> tfree = bm0, complex = 0)

    # ---- Observations ----
    # Total bevacizumab in mg/L: 1 nM = MW / 1e6 mg/L (MW in g/mol).
    Cc <- ctot * mw_bev / 1e6
    # Free VEGF165 in ng/L: 1 nM = MW (g/mol) ng/L.
    freeVegf <- tfree * mw_vegf165

    Cc       ~ prop(propSd)
    freeVegf ~ prop(propSd_freeVegf)
  })
}
