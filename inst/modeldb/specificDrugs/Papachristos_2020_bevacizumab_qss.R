Papachristos_2020_bevacizumab_qss <- function() {
  description <- "Quasi-steady-state target-mediated drug-disposition (TMDD QSS) model for IV bevacizumab and free VEGF-A in adults with metastatic colorectal cancer, with allometric weight scaling and ICAM-1 / VEGF-A genotype covariates (Papachristos 2020, Table 2)"
  reference <- "Papachristos A, Karatza E, Kalofonos H, Karalis V. Pharmacogenetics in Model-Based Optimization of Bevacizumab Therapy for Metastatic Colorectal Cancer. Int J Mol Sci. 2020;21(11):3753. doi:10.3390/ijms21113753"
  vignette <- "Papachristos_2020_bevacizumab_qss"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L (bevacizumab) and ng/L (free VEGF-A)")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed in source analysis; applied as power-form covariate on CL with reference weight 70 kg.",
      source_name        = "weight"
    ),
    SNP_ICAM1_RS1799969 = list(
      description        = "ICAM-1 rs1799969 mutant-allele-presence indicator (1 = heterozygous or homozygous mutant; 0 = homozygous wild-type)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (homozygous wild-type)",
      notes              = "Time-fixed per subject. Multiplicative effect on CL.",
      source_name        = "cat (Papachristos 2020 Table 2 narrative)"
    ),
    SNP_VEGFA_RS699947 = list(
      description        = "VEGF-A rs699947 mutant-allele-presence indicator (1 = heterozygous or homozygous mutant; 0 = homozygous wild-type)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (homozygous wild-type)",
      notes              = "Time-fixed per subject. Multiplicative effect on the QSS dissociation constant K_ss (mutant carriers have higher in-vivo affinity) and on baseline free VEGF-A BM0 (mutant carriers have lower baseline).",
      source_name        = "cat (Papachristos 2020 Table 2 narrative)"
    )
  )

  population <- list(
    n_subjects        = 46,
    n_studies         = 1,
    age_range         = "IQR 53-72 years",
    age_median        = "63 years",
    weight_range      = "IQR 64.15-81.75 kg",
    weight_median     = "74.5 kg",
    sex_female_pct    = 39,
    race_ethnicity    = "Greek (single-centre study); race not formally reported.",
    disease_state     = "Adults with metastatic colorectal cancer (mCRC); ECOG performance status 0-2.",
    dose_range        = "5 mg/kg IV every 2 weeks (76% of patients) or 7.5 mg/kg IV every 3 weeks (24%) with chemotherapy.",
    regions           = "Greece",
    n_observations    = "156 bevacizumab + 169 free VEGF-A serum concentrations",
    co_medication     = "BEV-FOLFIRI, BEV-FOLFOX, BEV-CapIRI, or BEV-CapOX.",
    notes             = "Single-centre prospective observational study (Papachristos 2020 sec. 2.1, 4.1)."
  )

  ini({
    # Structural PK parameters (reference 70 kg adult, all SNPs wild-type)
    lcl   <- log(0.344);  label("Linear clearance of free bevacizumab for a 70 kg wild-type adult (CL, L/day)")  # Table 2 row CL
    lv1   <- log(5.83);   label("Central volume of distribution (V1, L)")                                          # Table 2 row V1pop
    lq    <- log(0.136);  label("Inter-compartmental clearance (Q, L/day)")                                        # Table 2 row Qpop
    lv2   <- log(3.17);   label("Peripheral volume of distribution (V2, L)")                                       # Table 2 row V2pop

    # Target turnover (free VEGF-A)
    lkout <- log(0.116);  label("First-order elimination rate constant of free VEGF-A (kout, 1/day)")              # Table 2 row Koutpop
    lbm0  <- log(0.0137); label("Baseline (pre-dose) free VEGF-A concentration for a wild-type adult (BM0, nM)")   # Table 2 row BM0pop (= 616.5 ng/L)

    # Binding (QSS dissociation constant; same units as ctot, ttot)
    lkss  <- log(135);    label("Quasi-steady-state dissociation constant for bevacizumab-VEGF-A binding (Kss, nM)")  # Table 2 row KSSpop

    # Allometric / SNP covariate effects on log-CL
    allo_cl <- 1.01;  label("Allometric exponent on CL for log(WT/70) (unitless)")                                 # Table 2 row "log(weight/70) on CL"
    e_icam1_rs1799969_cl <- -0.33; label("ICAM-1 rs1799969 mutant effect on log-CL (unitless)")                    # Table 2 row "ICAM-1 rs1799969 mutant on CL"

    # SNP effects on Kss and BM0
    e_vegfa_rs699947_kss <-  1.22;  label("VEGF-A rs699947 mutant effect on log-Kss (unitless)")                   # Table 2 row "VEGF-A rs699947 mutant on KSS"
    e_vegfa_rs699947_bm0 <- -0.851; label("VEGF-A rs699947 mutant effect on log-BM0 (unitless)")                   # Table 2 row "VEGF-A rs699947 mutant on BM0"

    # IIV (exponential model; correlation rho(CL,Q) = -0.999 -> cov = -0.999 * 0.309 * 0.201 = -0.06205)
    etalcl + etalq ~ c(0.09548,
                       -0.06205, 0.04040)  # Table 2 omega_CL = 0.309, omega_Q = 0.201, p(Q,CL) = -0.999
    etalv1  ~ 0.02856  # Table 2 omega_V1 = 0.169
    etalv2  ~ 0.30803  # Table 2 omega_V2 = 0.555
    etalbm0 ~ 0.05760  # Table 2 omega_BM0 = 0.24

    # Residual error (proportional, two outputs)
    CcpropSd <- 0.253; label("Proportional residual error for bevacizumab (fraction)")  # Table 2 row sigma_BEVA
    CvpropSd <- 0.290; label("Proportional residual error for free VEGF-A (fraction)")  # Table 2 row sigma_VEGF
  })
  model({
    # ---- Unit-conversion constants ----
    # MW_BEV: bevacizumab molecular weight (g/mol). Standard reference value for the marketed
    # antibody; not reported numerically in Papachristos 2020 but required to convert mg dose
    # into the molar concentrations used by the binding equations. Documented as a model
    # assumption in the vignette's Assumptions and deviations section.
    mw_bev <- 149000

    # MW_VEGF (g/mol): back-calculated from Papachristos 2020 Table 2, which reports
    # BM0 = 0.0137 nM = 616.5 ng/L. MW_VEGF = 616.5 / 0.0137 = 45000 g/mol, consistent with
    # the homodimeric VEGF-A165 form (~22.5 kDa monomer).
    mw_vegf <- 45000

    # mg-to-nmol conversion factor for bevacizumab (1 mg = 1e6/MW nmol).
    nmol_per_mg_bev <- 1e6 / mw_bev

    # ---- Individual parameters (paper formulas, sec. 2.3) ----
    cl   <- exp(lcl + allo_cl * log(WT / 70) + e_icam1_rs1799969_cl * SNP_ICAM1_RS1799969 + etalcl)
    v1   <- exp(lv1 + etalv1)
    q    <- exp(lq + etalq)
    v2   <- exp(lv2 + etalv2)
    kout <- exp(lkout)
    bm0  <- exp(lbm0 + e_vegfa_rs699947_bm0 * SNP_VEGFA_RS699947 + etalbm0)
    kss  <- exp(lkss + e_vegfa_rs699947_kss * SNP_VEGFA_RS699947)

    # Zero-order target synthesis (kin = kout * BM0 at steady state with no drug)
    kin  <- kout * bm0

    # Complex elimination rate equals free-drug clearance per paper section 2.3
    # ("the elimination clearance of the bevacizumab-VEGF-A complex... was set equal to
    # the CL of free bevacizumab"), i.e. kint * V1 = CL, so kint = CL / V1.
    kint <- cl / v1

    # ---- QSS algebra (Gibiansky et al. 2008) ----
    # central tracks TOTAL bevacizumab amount (free + complex) in nmol;
    # total_target tracks TOTAL free+bound VEGF-A concentration in V1 (nM).
    ctot <- central / v1
    ttot <- total_target

    # Numerically stable discriminant: (Ctot + Ttot + Kss)^2 - 4 Ctot Ttot
    #                                 = (Ctot - Ttot)^2 + 2 (Ctot + Ttot) Kss + Kss^2
    qss_disc <- (ctot - ttot)^2 + 2 * (ctot + ttot) * kss + kss^2
    complex  <- ((ctot + ttot + kss) - sqrt(qss_disc)) / 2
    cfree    <- ctot - complex
    tfree    <- ttot - complex

    # ---- Dosing: convert mg input to nmol via bioavailability multiplier ----
    f(central) <- nmol_per_mg_bev

    # ---- ODEs ----
    # Central drug (free + complex; complex eliminated at rate kint * V1; peripheral has no binding)
    d/dt(central)     <- -cl * cfree - q * cfree + q * (peripheral1 / v2) - kint * complex * v1
    d/dt(peripheral1) <-               q * cfree - q * (peripheral1 / v2)
    # Total target (free + complex; only present in central per QSS-TMDD convention)
    d/dt(total_target) <- kin - kout * tfree - kint * complex
    total_target(0)    <- bm0  # steady-state baseline (no drug -> tfree = bm0, complex = 0)

    # ---- Observations ----
    # Bevacizumab observation: total drug concentration in mg/L
    # nM -> mg/L conversion: c_nM * MW_BEV (g/mol) * 1e-9 (g/ng/L conversion not needed here)
    #   1 nM = 1e-9 mol/L * MW (g/mol) = MW * 1e-9 g/L = MW / 1e6 mg/L
    Cc <- ctot * mw_bev / 1e6
    Cc ~ prop(CcpropSd)

    # Free VEGF-A observation: ng/L (paper unit)
    # 1 nM = MW (g/mol) * 1e-9 g/L = MW ng/L
    Cv <- tfree * mw_vegf
    Cv ~ prop(CvpropSd)
  })
}
