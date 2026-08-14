Michelet_2025_BI754111_mpbpk <- function() {
  description <- paste(
    "QSP / minimal-PBPK. Two-pore minimal physiologically-based pharmacokinetic",
    "platform for the anti-LAG-3 IgG monoclonal antibody BI 754111, developed to",
    "predict intratumor exposure and LAG-3 receptor occupancy from 89Zr-immuno-PET",
    "biodistribution data in patients with NSCLC or HNSCC. Couples a two-compartment",
    "plasma popPK (parallel linear and saturable clearance) to a Shah-Betts (2012)",
    "physiologic tumor tissue compartment (vascular, endosomal and interstitial",
    "sub-spaces with FcRn recycling) plus a draining lymph-node space, with dynamic",
    "BI 754111-LAG-3 binding in both blood and tumor, internalization of the complex",
    "into the T-cellular space, removal of internalized complex by T-cell turnover,",
    "and an intracellular LAG-3 pool that recycles receptor back to the cell surface.",
    "The 89Zr-labeled tracer and the unlabeled drug are carried as two parallel",
    "species (labeled states use the _zr89 suffix) because residualizing 89Zr stays",
    "inside the T cell after internalization and keeps contributing to the PET signal",
    "whereas unlabeled antibody is lost from the system. Extravasation from the tumor",
    "vascular into the interstitial space follows the Li-Shah two-pore formalism",
    "rather than a single vascular reflection coefficient. Parameters are the final",
    "average two-pore model with internal LAG-3 pool (Michelet 2025 Table 2)."
  )
  reference <- paste(
    "Michelet R, Petersson K, Huisman MC, Menke-van der Houven van Oordt CW,",
    "Miedema IHC, Thiele A, Montaseri G, Perez-Pitarch A, Busse D. A minimal",
    "physiologically-based pharmacokinetic modeling platform to predict intratumor",
    "exposure and receptor occupancy of an anti-LAG-3 monoclonal antibody.",
    "CPT Pharmacometrics Syst Pharmacol. 2025;14(3):460-473. doi:10.1002/psp4.13285.",
    "Model code in Appendix S1 (mrgsolve), key equations in Appendix S3, parameter",
    "table in Appendix S4 Table S1. The mPBPK backbone is adapted from Lindauer A et",
    "al. CPT Pharmacometrics Syst Pharmacol. 2017;6(1):11-20 doi:10.1002/psp4.12130",
    "(see modellib('Lindauer_2017_pembrolizumab')); the tumor tissue structure and",
    "physiologic constants follow Shah DK, Betts AM. J Pharmacokinet Pharmacodyn.",
    "2012;39:67-86. The two-pore formalism follows Li Z, Shah DK. J Pharmacokinet",
    "Pharmacodyn. 2019;46:305-318. The plasma popPK and the 89Zr-immuno-PET",
    "biodistribution data are from Miedema IHC et al. Eur J Nucl Med Mol Imaging.",
    "2023;50:2068-2080 (trial NCT03780725).",
    sep = " "
  )
  vignette <- "Michelet_2025_BI754111_mpbpk"

  # Both species are dosed: the unlabeled therapeutic mass dose into `central`
  # and the 89Zr-labeled tracer dose into `central_zr89`. The latter is outside
  # the depot/central heuristic, so it is declared explicitly.
  dosing <- c("central", "central_zr89")

  # Every state below other than `central` / `peripheral1` / `lnode` is a
  # paper-mechanistic state of the Shah-Betts / Lindauer tumor structure or of
  # the 89Zr-labeled tracer arm. The `_zr89` suffix marks the radiolabeled
  # species; it is deliberately declared paper-specific rather than registered
  # as a shared metabolite/sibling-drug suffix (see the vignette's Assumptions
  # and deviations section).
  paper_specific_compartments <- c(
    "tumor_vs", "tumor_is", "tumor_es_ub", "tumor_es_b",
    "complex_tumor", "complex_tumor_int", "complex_blood", "complex_blood_int",
    "target_tumor", "target_tumor_int", "target_blood", "target_blood_int",
    "fcrn",
    "central_zr89", "peripheral1_zr89", "lnode_zr89",
    "tumor_vs_zr89", "tumor_is_zr89", "tumor_es_ub_zr89", "tumor_es_b_zr89",
    "complex_tumor_zr89", "complex_tumor_int_zr89",
    "complex_blood_zr89", "complex_blood_int_zr89"
  )

  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Note that the Michelet 2025 Appendix S1 implementation
  # deliberately mixes amount states (nmol) and concentration states (nM); that
  # mixture is preserved here so the equations match the published code.
  compartmentData <- list(
    central                = list(analyte = "BI 754111", units = "nmol", specimen = "plasma", verified = TRUE),
    peripheral1            = list(analyte = "BI 754111", units = "nmol", specimen = "plasma", verified = TRUE),
    lnode                  = list(analyte = "BI 754111", units = "nmol", specimen = "lymph", verified = TRUE),
    tumor_vs               = list(analyte = "BI 754111", units = "nmol", specimen = "tumor", verified = TRUE),
    tumor_es_ub            = list(analyte = "BI 754111", units = "nM", specimen = "endosome", verified = TRUE),
    tumor_es_b             = list(analyte = "BI 754111", units = "nM", specimen = "endosome", verified = TRUE),
    tumor_is               = list(analyte = "BI 754111", units = "nmol", specimen = "tumor", verified = TRUE),
    complex_tumor          = list(analyte = "BI 754111-LAG-3 complex", units = "nM", specimen = "tumor", verified = TRUE),
    complex_tumor_int      = list(analyte = "BI 754111-LAG-3 complex", units = "nmol", specimen = "tumor", verified = TRUE),
    complex_blood          = list(analyte = "BI 754111-LAG-3 complex", units = "nM", specimen = "whole blood", verified = TRUE),
    complex_blood_int      = list(analyte = "BI 754111-LAG-3 complex", units = "nmol", specimen = "whole blood", verified = TRUE),
    target_tumor           = list(analyte = "LAG-3", units = "nmol", specimen = "tumor", verified = TRUE),
    target_tumor_int       = list(analyte = "LAG-3", units = "nM", specimen = "tumor", verified = TRUE),
    target_blood           = list(analyte = "LAG-3", units = "nmol", specimen = "whole blood", verified = TRUE),
    target_blood_int       = list(analyte = "LAG-3", units = "nM", specimen = "whole blood", verified = TRUE),
    fcrn                   = list(analyte = "FcRn", units = "nM", specimen = "endosome", verified = TRUE),
    central_zr89           = list(analyte = "89Zr-labeled BI 754111", units = "nmol", specimen = "plasma", verified = TRUE),
    peripheral1_zr89       = list(analyte = "89Zr-labeled BI 754111", units = "nmol", specimen = "plasma", verified = TRUE),
    lnode_zr89             = list(analyte = "89Zr-labeled BI 754111", units = "nmol", specimen = "lymph", verified = TRUE),
    tumor_vs_zr89          = list(analyte = "89Zr-labeled BI 754111", units = "nmol", specimen = "tumor", verified = TRUE),
    tumor_es_ub_zr89       = list(analyte = "89Zr-labeled BI 754111", units = "nM", specimen = "endosome", verified = TRUE),
    tumor_es_b_zr89        = list(analyte = "89Zr-labeled BI 754111", units = "nM", specimen = "endosome", verified = TRUE),
    tumor_is_zr89          = list(analyte = "89Zr-labeled BI 754111", units = "nmol", specimen = "tumor", verified = TRUE),
    complex_tumor_zr89     = list(analyte = "89Zr-BI 754111-LAG-3 complex", units = "nM", specimen = "tumor", verified = TRUE),
    complex_tumor_int_zr89 = list(analyte = "89Zr-BI 754111-LAG-3 complex", units = "nmol", specimen = "tumor", verified = TRUE),
    complex_blood_zr89     = list(analyte = "89Zr-BI 754111-LAG-3 complex", units = "nM", specimen = "whole blood", verified = TRUE),
    complex_blood_int_zr89 = list(analyte = "89Zr-BI 754111-LAG-3 complex", units = "nmol", specimen = "whole blood", verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species         = "human",
    n_subjects      = 6L,
    n_studies       = 1L,
    age_range       = "not reported in Michelet 2025 or its appendices",
    weight_range    = "not reported in Michelet 2025 or its appendices",
    sex_female_pct  = NA_real_,
    race_ethnicity  = "not reported",
    disease_state   = paste(
      "Advanced solid tumors having progressed on anti-PD-1 based treatment:",
      "non-small cell lung cancer (n = 4) and head and neck squamous cell",
      "carcinoma (n = 2). All participants continued anti-PD-1 treatment with",
      "ezabenlimab during the imaging study."
    ),
    dose_range      = paste(
      "4 mg 89Zr-labeled BI 754111 tracer dose in cycle 1; two weeks later a",
      "40 mg or 600 mg unlabeled BI 754111 mass dose followed by a second 4 mg",
      "tracer dose. The underlying popPK was developed on 4-600 mg IV doses in",
      "49 patients (phase I NCT03156114)."
    ),
    regions         = "Netherlands (Amsterdam UMC); trial NCT03780725",
    notes           = paste(
      "PET scans were acquired < 2, 90 +/- 1 and 138 +/- 1 h after the cycle-1",
      "tracer injection and at 90 +/- 1 and 138 +/- 1 h after the cycle-2 tracer",
      "injection. The mPBPK layer is a deterministic platform model: it was",
      "calibrated by minimizing a weighted combination of average fold error,",
      "absolute average fold error and residual sum of squared errors (Appendix",
      "S2 Equations 3-6) rather than by nonlinear mixed-effects estimation, so no",
      "interindividual variability and no residual-error model are reported for",
      "the mPBPK parameters. The authors state explicitly that quantifying",
      "between-patient variability with a nonlinear mixed-effects approach is",
      "future work. Plasma PK is driven by empirical Bayes estimates of the",
      "Miedema 2023 popPK model; Appendix S4 Table S1 reports the six individual",
      "CL and V1 estimates but no population typical value for those two",
      "parameters, so the defaults in this file are the medians of the six",
      "published individual estimates (see the validation vignette). The",
      "structural parameters that were re-estimated for the final model",
      "(T_multi, N_Tcell, N_LAG-3,TC, K_deg,TC, K_rec,LAG-3, K_int,LAG-3) are the",
      "average of the eight best-performing model runs reported in Table 2."
    )
  )

  ini({
    # ======================================================================
    # Drug and physiologic constants
    # ======================================================================
    mw_mab   <- fixed(149000); label("BI 754111 molecular weight (g/mol)")                                            # Table S1 row MW = 149,000 g/mol (Boehringer Ingelheim data on file)
    v_blood  <- fixed(5);      label("Blood volume (L; = 5,000,000 uL)")                                              # Table S1 row V_blood = 5,000,000 uL
    v_lnode  <- fixed(0.274);  label("Volume of the (benign) draining lymph-node compartment (L; = 274,000 uL = 274 mL)")  # Table S1 row V_LN = 274,000 uL. Appendix S1 declares V_LN = 274 "mL" then divides by 1e6, giving 2.74e-4 L; its own iL_LN comment (PLQ*V_LN*1e-3 = 3.4798 L/h) and the Shah-Betts human lymph-node volume both confirm 0.274 L. See vignette Errata.

    # ======================================================================
    # Plasma popPK (Miedema 2023, reported in Michelet 2025 Table S1).
    # Two-compartment with parallel linear and saturable clearance.
    # Table S1 reports population values for V2, Q, CL_SAT and C50, but only
    # the six individual empirical Bayes estimates for CL and V1; the defaults
    # below are the medians of those six values (footnotes a and c):
    #   CL (L/h): 0.0129, 0.0305, 0.0131, 0.00175, 0.0104, 0.0113 -> median 0.0121
    #   V1 (L):   2.64, 3.25, 2.50, 2.74, 2.47, 3.54              -> median 2.69
    # All plasma-PK parameters are inherited from an upstream popPK fit and are
    # therefore fixed here.
    # ======================================================================
    lcl     <- fixed(log(0.0121));  label("Linear clearance CL (L/h)")                                                # Table S1 footnote c: median of the six individual EBEs (0.00175-0.0305 L/h)
    lvc     <- fixed(log(2.69));    label("Central (blood) volume V1 (L)")                                            # Table S1 footnote a: median of the six individual EBEs (2.47-3.54 L)
    lq      <- fixed(log(0.0353));  label("Inter-compartmental clearance Q (L/h)")                                    # Table S1 row Q = 0.0353 L/h (population estimate)
    lvp     <- fixed(log(2.25));    label("Peripheral volume V2 (L)")                                                 # Table S1 row V2 = 2.25 L (population estimate)
    lclsat  <- fixed(log(0.0489));  label("Saturable (non-linear) clearance CL_SAT (L/h)")                            # Table S1 row CL_SAT = 0.0489 L/h (population estimate)
    lc50    <- fixed(log(14.497));  label("Concentration giving 50% of the saturable elimination C50 (nM; = 2160 ug/L / 149000 g/mol)")  # Table S1 row C50 = 2160 ug/L, converted to nM with MW

    # ======================================================================
    # Tumor tissue structure (Shah & Betts 2012, via Table S1). Identical to
    # the values used by modellib('Lindauer_2017_pembrolizumab').
    # ======================================================================
    f_v_es     <- fixed(0.005);  label("Endosomal-space volume as fraction of total tumor volume (unitless)")         # Table S1 row V_es,per = 0.50% of TV
    f_v_is     <- fixed(0.55);   label("Interstitial-space volume as fraction of total tumor volume (unitless)")      # Table S1 row V_is,per = 55% of TV
    f_v_vs     <- fixed(0.07);   label("Vascular-space volume as fraction of total tumor volume (unitless)")          # Table S1 row V_vs,per = 7% of TV
    plq_norm   <- fixed(12.7);   label("Tumor plasma flow per unit tissue volume (L/h/L)")                            # Table S1 row PLQ_rel = 12.7 L/h/L
    f_lymph    <- fixed(0.002);  label("Lymph flow as fraction of tumor plasma flow (unitless; = 0.20% of PLQ)")      # Table S1 row L_per = 0.20% of PLQ
    clup_norm  <- fixed(0.0366); label("Endosomal pinocytosis per unit endosomal-space volume (L/h/L)")               # Appendix S1 param CLup_rel = 0.0366 L/h/L (Table S1 quotes the rounded 0.04)
    kdeg_endo  <- fixed(42.9);   label("Endosomal degradation rate constant of free antibody (1/h)")                  # Table S1 row K_deg = 42.9 1/h
    v_ref_is   <- fixed(0.2);    label("Lymph / interstitial reflection coefficient (unitless)")                      # Table S1 row V_ref,is = 0.2
    fcrn_init  <- fixed(49800);  label("Initial endosomal FcRn concentration (nM; = 49.8 uM)")                        # Table S1 row FcRni = 49.8 uM
    fr_recycle <- fixed(0.715);  label("Fraction of FcRn-bound antibody recycled to the vascular space (unitless)")   # Table S1 row FR = 0.715
    kon_fcrn   <- fixed(0.792);  label("FcRn-antibody association rate constant (1/(nM*h); = 792 x 1e6 /M/h)")        # Table S1 row K_on,FcRn = 792 x 1e6 /M/h -> x1e6/1e9
    koff_fcrn  <- fixed(23.9);   label("FcRn-antibody dissociation rate constant (1/h)")                              # Table S1 row K_off,FcRn = 23.9 1/h

    # ======================================================================
    # Two-pore extravasation (Li & Shah 2019, implemented in Appendix S1).
    # These replace the single vascular reflection coefficient V_ref of the
    # starting model; the resulting effective V_ref for BI 754111 is ~0.68
    # (Michelet 2025 Results, "Inclusion of two-pore formalism").
    # ======================================================================
    r_pore_l    <- fixed(22.85); label("Large-pore radius (nm)")                                                      # Appendix S1 param r_L = 22.85 nm (Li & Shah 2019)
    r_pore_s    <- fixed(4.44);  label("Small-pore radius (nm)")                                                      # Appendix S1 param r_S = 4.44 nm (Li & Shah 2019)
    alpha_l     <- fixed(0.042); label("Fractional hydraulic conductance of large pores (unitless)")                  # Appendix S1 param alpha_L = 0.042
    alpha_s     <- fixed(0.958); label("Fractional hydraulic conductance of small pores (unitless)")                  # Appendix S1 param alpha_S = 0.958
    sigma_alb_l <- fixed(0.108); label("Albumin reflection coefficient for large pores (unitless)")                   # Appendix S1 param sigma_alb_L = 0.108
    sigma_alb_s <- fixed(0.95);  label("Albumin reflection coefficient for small pores (unitless)")                   # Appendix S1 param sigma_alb_S = 0.95
    delta_p     <- fixed(10);    label("Osmotic pressure across the tissue vascular endothelial membrane (mmHg)")     # Appendix S1 param delta_P = 10 mmHg
    st_force    <- fixed(1);     label("Starling force (mmHg)")                                                       # Appendix S1 param St_force = 1 mmHg
    rgas        <- fixed(62.363); label("Gas constant (L*mmHg/K/mol)")                                                # Appendix S1 param Rgas = 62.363
    temp_k      <- fixed(310);   label("Body temperature (K)")                                                        # Appendix S1 param Temp = 310 K

    # ======================================================================
    # BI 754111-LAG-3 binding (Boehringer Ingelheim in vitro, via Table S1).
    # ======================================================================
    kon_lag3  <- fixed(5.2);     label("BI 754111-LAG-3 association rate constant (1/(nM*h); = 5200 x 1e6 /M/h)")     # Table S1 row K_on,LAG-3 = 5200 x 1e6 /M/h -> x1e6/1e9
    koff_lag3 <- fixed(0.144);   label("BI 754111-LAG-3 dissociation rate constant (1/h)")                            # Table S1 row K_off,LAG-3 = 0.144 1/h
    kout_lag3 <- fixed(0.00246); label("LAG-3 turnover (degradation) rate constant (1/h); K_out = K_deg,LAG-3")       # Table S1 row K_deg,LAG-3 = 0.00246 1/h; Appendix S3 states K_out was assumed equal to K_deg,LAG-3

    # ======================================================================
    # Final model parameters -- the average of the eight best-performing runs
    # of the two-pore model with internal LAG-3 pool (Michelet 2025 Table 2,
    # column "Final average two-pore model (with LAG-3-pool)"). These are the
    # parameters that were re-estimated against the PET biodistribution data;
    # the starting-model values are given in the comment on each line.
    # ======================================================================
    n_tcell    <- fixed(6430);    label("Activated (CD4+/CD8+) T cells per uL of blood")                              # Table 2 final average two-pore model N_Tcell = 6430 (starting model 2000)
    n_lag3_tc  <- fixed(1330);    label("LAG-3 receptors per activated T cell")                                       # Table 2 final average two-pore model N_LAG-3,TC = 1330 (starting model 511)
    tmulti     <- fixed(16.2);    label("Initial ratio of LAG-3 concentration in tumor versus blood (unitless)")      # Table 2 final average two-pore model T_multi = 16.2 (starting model 4.3)
    kint_lag3  <- fixed(0.466);   label("Internalization rate of the BI 754111-LAG-3 complex into T cells (1/h)")     # Table 2 final average two-pore model K_int,LAG-3 = 0.466 1/h (starting model 0.276)
    krec_lag3  <- fixed(0.453);   label("Fraction of intracellular LAG-3 recycled to the cell surface per hour")      # Table 2 final average two-pore model K_rec,LAG-3 = 0.453 (starting model 0.30)
    kdeg_tc    <- fixed(0.00725); label("Removal of internalized complex and receptor by T-cell turnover (1/h)")      # Table 2 final average two-pore model K_deg,TC = 0.00725 1/h (starting model 0.00963)

    # ======================================================================
    # Tumor size. The tumor volume is derived from the baseline sum of longest
    # diameters via the ellipsoid shape factor coded in Appendix S1; it is
    # constant over the simulation (dTV/dt = 0 -- no tumor-growth sub-model).
    # ======================================================================
    lw0_sld <- log(7.5); label("Baseline tumor size as sum of longest diameters (cm)")                                # Appendix S1 param W0_SLD = 7.5 cm

    # ======================================================================
    # Residual error. The mPBPK layer is deterministic and no residual-error
    # model is reported anywhere in Michelet 2025 or its appendices; a nominal
    # proportional term is carried so the model has a defined endpoint.
    # ======================================================================
    propSd <- fixed(0.001); label("Residual proportional error (nominal; deterministic platform model)")              # Not reported in the source; nominal placeholder, see vignette Errata
  })

  model({
    # ==================================================================
    # Internal unit conventions (as in Michelet 2025 Appendix S1)
    #   Time                     h
    #   Volumes                  L
    #   Antibody amounts         nmol
    #   Antibody concentrations  nM (= nmol/L)
    #   All Kon                  1/(nM*h)
    #   All Koff, Kdeg, Kint     1/h
    # Note that Appendix S1 carries the endosomal, complex and internal-pool
    # states as CONCENTRATIONS while the plasma, lymph, vascular and
    # interstitial states are AMOUNTS. That mixture is reproduced exactly so
    # the equations match the published code.
    # ==================================================================

    cl    <- exp(lcl)
    vc    <- exp(lvc)
    q     <- exp(lq)
    vp    <- exp(lvp)
    clsat <- exp(lclsat)
    c50   <- exp(lc50)
    w0sld <- exp(lw0_sld)

    # ---- micro-constants of the plasma popPK ----
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    km  <- clsat / vc

    # ---- unit conversion: a dose in mg lands as nmol ----
    mg_to_nmol <- 1e6 / mw_mab
    f(central)      <- mg_to_nmol
    f(central_zr89) <- mg_to_nmol

    # ==================================================================
    # Tumor geometry. Appendix S1 converts the baseline sum of longest
    # diameters to a volume with an ellipsoid shape factor built from
    # 1:1 and 1:5 axis ratios, then converts cm^3 -> mm^3 (= uL).
    # ==================================================================
    shape   <- (((1 / 5)^2)^(1 / 3) + ((1 / 1)^2)^(1 / 3)) / 2
    gam     <- 3.141593 / 6 * shape^3
    w0      <- gam * w0sld^3 * 1000
    v_total <- w0 / 1e6
    v_vs    <- v_total * f_v_vs
    v_is    <- v_total * f_v_is
    v_es    <- v_total * f_v_es

    # ---- physiologic flows in the tumor and the lymph node (L/h) ----
    plq    <- plq_norm * v_total
    lymph  <- plq * f_lymph
    clup   <- clup_norm * v_es
    plq_ln <- plq_norm * v_lnode
    l_ln   <- plq_ln * f_lymph

    # ==================================================================
    # Two-pore extravasation (Li & Shah 2019; Appendix S1 [ode] block).
    # a_e is the hydrodynamic radius of the antibody; sigma_S / sigma_L the
    # per-pore reflection coefficients; A_Ao_* the fractional accessible
    # pore areas; XJ the isogravimetric flow; Pe_* the Peclet numbers.
    # The resulting CL_vs_is replaces (1 - v_ref) * lymph of the one-pore
    # model and corresponds to an effective v_ref of about 0.68.
    # ==================================================================
    avogadro <- 6.0221415e23
    a_e      <- 0.0483 * mw_mab^0.386
    sigma_s  <- 1 - 0.8489 * exp(-0.00004 * mw_mab)
    sigma_l  <- 0.000035 * mw_mab^0.717
    a_ao_s   <- 0.2352 * exp(-0.00008295 * mw_mab) + 0.7767 * exp(-0.00053095 * mw_mab)
    a_ao_l   <- 0.3429 * exp(-0.00012175 * mw_mab) + 0.6571 * exp(-0.00000421 * mw_mab)
    xp       <- (rgas * temp_k / (6 * 3.141593 * avogadro)) * (8 / st_force) * 1e24
    xp_s     <- xp * a_ao_s * alpha_s / (a_e * r_pore_s^2)
    xp_l     <- xp * a_ao_l * alpha_l / (a_e * r_pore_l^2)
    xj       <- alpha_l * alpha_s * (sigma_alb_s - sigma_alb_l) * delta_p / st_force
    pe_s     <- (-xj + alpha_s) * (1 - sigma_s) / xp_s
    pe_l     <- (xj + alpha_l) * (1 - sigma_l) / xp_l

    # ---- concentrations (nM) of the unlabeled and labeled species ----
    Cp    <- central     / vc
    Cvs   <- tumor_vs    / v_vs
    Cis   <- tumor_is    / v_is
    Clym  <- lnode       / v_lnode
    Cp_h  <- central_zr89  / vc
    Cvs_h <- tumor_vs_zr89 / v_vs
    Cis_h <- tumor_is_zr89 / v_is
    Clym_h <- lnode_zr89   / v_lnode

    # ---- two-pore clearance from vascular into interstitial space (L/h) ----
    cis_cvs_ratio <- (Cis + Cis_h + 1e-7) / (Cvs + Cvs_h + 1e-7)
    cl_conv <- (-xj + alpha_s) * (1 - sigma_s) * lymph +
      (xj + alpha_l) * (1 - sigma_l) * lymph
    cl_diff <- xp_s * lymph * (1 - cis_cvs_ratio) * pe_s / (exp(pe_s) - 1) +
      xp_l * lymph * (1 - cis_cvs_ratio) * pe_l / (exp(pe_l) - 1)
    cl_vs_is <- cl_conv + cl_diff

    # ==================================================================
    # LAG-3 receptor baselines. The number of receptors in blood follows
    # from the T-cell count and the receptors per T cell; the tumor amount
    # is T_multi times the blood concentration in the interstitial volume.
    # ==================================================================
    m_lag3_bi <- n_lag3_tc * n_tcell * (v_blood * 1e6) / avogadro * 1e9
    c_lag3_bi <- m_lag3_bi / v_blood
    v_is_in   <- f_v_is * w0 / 1e6
    m_lag3_ti <- c_lag3_bi * tmulti * v_is_in
    kin_t     <- m_lag3_ti * kout_lag3
    kin_b     <- m_lag3_bi * kout_lag3

    # ---- recycling rate constant from the in vitro percentage per hour ----
    krec <- log(2) / (0.5 / krec_lag3)

    # ---- free (unoccupied) receptor concentrations (nM) ----
    C_lag3_t <- target_tumor / v_is
    C_lag3_b <- target_blood / v_blood
    free_t   <- C_lag3_t - complex_tumor - complex_tumor_zr89
    free_b   <- C_lag3_b - complex_blood - complex_blood_zr89

    # ==================================================================
    # 89Zr-labeled BI 754111 ("hot"). Identical structure to the unlabeled
    # arm except that internalized complex residualizes inside the T cell
    # and keeps contributing to the tumor signal.
    # ==================================================================
    d/dt(central_zr89) <-
      -kel * central_zr89 + k21 * peripheral1_zr89 - k12 * central_zr89 -
      km * central_zr89 / (1 + (Cp_h + Cp) / c50) -
      plq * Cp_h + (plq - lymph) * Cvs_h + l_ln * Clym_h

    d/dt(peripheral1_zr89) <- k12 * central_zr89 - k21 * peripheral1_zr89

    d/dt(lnode_zr89) <- (1 - v_ref_is) * lymph * Cis_h - l_ln * Clym_h

    d/dt(tumor_vs_zr89) <-
      plq * Cp_h - (plq - lymph) * Cvs_h - cl_vs_is * Cvs_h -
      clup * Cvs_h + clup * fr_recycle * tumor_es_b_zr89

    d/dt(tumor_es_ub_zr89) <-
      clup / v_es * (Cvs_h + Cis_h) - kon_fcrn * tumor_es_ub_zr89 * fcrn +
      koff_fcrn * tumor_es_b_zr89 - kdeg_endo * tumor_es_ub_zr89

    d/dt(tumor_es_b_zr89) <-
      kon_fcrn * tumor_es_ub_zr89 * fcrn - koff_fcrn * tumor_es_b_zr89 -
      clup / v_es * tumor_es_b_zr89

    d/dt(tumor_is_zr89) <-
      cl_vs_is * Cvs_h - (1 - v_ref_is) * lymph * Cis_h - clup * Cis_h +
      clup * (1 - fr_recycle) * tumor_es_b_zr89 -
      kon_lag3 * Cis_h * v_is * free_t + koff_lag3 * complex_tumor_zr89 * v_is

    d/dt(complex_tumor_zr89) <-
      kon_lag3 * Cis_h * free_t - koff_lag3 * complex_tumor_zr89 -
      kint_lag3 * complex_tumor_zr89

    # Appendix S3 Equation 12
    d/dt(complex_tumor_int_zr89) <-
      kint_lag3 * complex_tumor_zr89 * v_is - kdeg_tc * complex_tumor_int_zr89

    d/dt(complex_blood_zr89) <-
      kon_lag3 * Cp_h * free_b - koff_lag3 * complex_blood_zr89 -
      kint_lag3 * complex_blood_zr89

    d/dt(complex_blood_int_zr89) <-
      kint_lag3 * complex_blood_zr89 * v_blood - kdeg_tc * complex_blood_int_zr89

    # ==================================================================
    # Unlabeled BI 754111 ("cold").
    # ==================================================================
    d/dt(central) <-
      -kel * central + k21 * peripheral1 - k12 * central -
      km * central / (1 + (Cp_h + Cp) / c50) -
      plq * Cp + (plq - lymph) * Cvs + l_ln * Clym

    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    d/dt(lnode) <- (1 - v_ref_is) * lymph * Cis - l_ln * Clym

    d/dt(tumor_vs) <-
      plq * Cp - (plq - lymph) * Cvs - cl_vs_is * Cvs -
      clup * Cvs + clup * fr_recycle * tumor_es_b

    d/dt(tumor_es_ub) <-
      clup / v_es * (Cvs + Cis) - kon_fcrn * tumor_es_ub * fcrn +
      koff_fcrn * tumor_es_b - kdeg_endo * tumor_es_ub

    d/dt(tumor_es_b) <-
      kon_fcrn * tumor_es_ub * fcrn - koff_fcrn * tumor_es_b -
      clup / v_es * tumor_es_b

    d/dt(tumor_is) <-
      cl_vs_is * Cvs - (1 - v_ref_is) * lymph * Cis - clup * Cis +
      clup * (1 - fr_recycle) * tumor_es_b -
      kon_lag3 * Cis * v_is * free_t + koff_lag3 * complex_tumor * v_is

    d/dt(complex_tumor) <-
      kon_lag3 * Cis * free_t - koff_lag3 * complex_tumor -
      kint_lag3 * complex_tumor

    d/dt(complex_tumor_int) <-
      kint_lag3 * complex_tumor * v_is - kdeg_tc * complex_tumor_int

    d/dt(complex_blood) <-
      kon_lag3 * Cp * free_b - koff_lag3 * complex_blood -
      kint_lag3 * complex_blood

    d/dt(complex_blood_int) <-
      kint_lag3 * complex_blood * v_blood - kdeg_tc * complex_blood_int

    # ==================================================================
    # LAG-3 receptor turnover with the intracellular pool.
    # Appendix S3 Equations 13 (pool) and 14 (surface receptor); the
    # internalization terms of Equation 14 both remove surface receptor,
    # matching the Appendix S1 code (the supplement prints the labeled and
    # unlabeled terms with opposite signs, which is a transcription error).
    # ==================================================================
    d/dt(target_tumor) <-
      kin_t - kout_lag3 * target_tumor -
      kint_lag3 * complex_tumor * v_is - kint_lag3 * complex_tumor_zr89 * v_is +
      krec * target_tumor_int * v_is

    d/dt(target_tumor_int) <-
      kint_lag3 * complex_tumor + kint_lag3 * complex_tumor_zr89 -
      kdeg_tc * target_tumor_int - krec * target_tumor_int

    d/dt(target_blood) <-
      kin_b - kout_lag3 * target_blood -
      kint_lag3 * complex_blood * v_blood - kint_lag3 * complex_blood_zr89 * v_blood +
      krec * target_blood_int * v_blood

    d/dt(target_blood_int) <-
      kint_lag3 * complex_blood + kint_lag3 * complex_blood_zr89 -
      kdeg_tc * target_blood_int - krec * target_blood_int

    # ---- endosomal FcRn shared by both species ----
    d/dt(fcrn) <-
      -kon_fcrn * tumor_es_ub_zr89 * fcrn + koff_fcrn * tumor_es_b_zr89 +
      clup / v_es * tumor_es_b_zr89 -
      kon_fcrn * tumor_es_ub * fcrn + koff_fcrn * tumor_es_b +
      clup / v_es * tumor_es_b

    # ---- initial conditions ----
    # Appendix S3 Equations 15-16: about half of the LAG-3 molecules are
    # retained intracellularly at baseline.
    target_tumor(0)     <- m_lag3_ti
    target_blood(0)     <- m_lag3_bi
    target_tumor_int(0) <- 0.5 * m_lag3_ti / v_is_in
    target_blood_int(0) <- 0.5 * c_lag3_bi
    fcrn(0)             <- fcrn_init

    # ==================================================================
    # Outputs. Concentrations are reported in ug/mL as in Appendix S1
    # [table]. The labeled tumor concentration includes the residualized
    # internalized complex (which the PET camera still sees); the unlabeled
    # tumor concentration does not, because unlabeled antibody is lost from
    # the system on internalization.
    # ==================================================================
    Cc <- (central / vc) * mw_mab / 1e6

    Cc_zr89 <- (central_zr89 / vc) * mw_mab / 1e6

    Ctumor <- ((tumor_vs + tumor_is + complex_tumor * v_is +
      tumor_es_b * v_es + tumor_es_ub * v_es) / v_total) * mw_mab / 1e6

    Ctumor_zr89 <- ((tumor_vs_zr89 + tumor_is_zr89 + complex_tumor_zr89 * v_is +
      complex_tumor_int_zr89 + tumor_es_b_zr89 * v_es +
      tumor_es_ub_zr89 * v_es) / v_total) * mw_mab / 1e6

    ROtumor <- 100 * (complex_tumor + complex_tumor_zr89) / C_lag3_t
    ROblood <- 100 * (complex_blood + complex_blood_zr89) / C_lag3_b

    Cc ~ prop(propSd)
  })
}
