Hu_2026_pemetrexed_osimertinib_mouse_qsp <- function() {
  description <- paste(
    "QSP. Preclinical (mouse, female BALB/c nude with subcutaneous HCC827",
    "EGFR-mutant NSCLC xenograft).",
    "Mechanistic QSP-PK-PD model for the sequence-dependent synergy between",
    "pemetrexed (PEM) and the third-generation EGFR-TKI osimertinib (OSI).",
    "Five coupled modules: two-compartment PEM PK (intraperitoneal),",
    "two-compartment OSI PK (oral), a folate module in which PEM accelerates",
    "degradation of the folate-metabolizing enzyme pool and thereby depletes",
    "folate, an EGFR module in which OSI depletes EGFR signal while damaged",
    "tumour cells drive a compensatory rebound of EGFR synthesis, and a",
    "Simeoni-type tumour-growth-inhibition module (one proliferating pool plus",
    "a three-state damaged-cell transit chain). Folate depletion accelerates",
    "the proliferating-to-damaged transition; EGFR suppression both slows",
    "proliferation and (via the Bim surrogate 1 - EGFR) accelerates apoptosis",
    "of damaged cells, while EGFR-driven G1 arrest antagonises the",
    "proliferating-to-damaged transition. Sequence dependence emerges from the",
    "model structure and dosing schedule alone -- no schedule-specific",
    "interaction parameter is fitted. An unperturbed (drug-free) twin tumour",
    "chain is carried alongside so real-time TGI% can be computed.",
    "Parameter values from Hu 2026 supplementary Tables S17-S21."
  )
  reference <- paste(
    "Hu K, Lin Y, Ji H, Yuan T, Xia Y, Yang J.",
    "A Mechanistic Pharmacokinetic/Pharmacodynamic Model for",
    "Sequence-Dependent Synergy in Pemetrexed-Osimertinib Combinations",
    "Against Non-Small Cell Lung Cancer (NSCLC): Translational Insights.",
    "Pharmaceutics. 2026;18(4):408.",
    "doi:10.3390/pharmaceutics18040408.",
    sep = " "
  )
  vignette <- "Hu_2026_pemetrexed_osimertinib_mouse_qsp"

  # Two separately dosed drugs: pemetrexed enters `depot` (intraperitoneal) and
  # osimertinib enters `depot_osimertinib` (oral gavage). Declared explicitly
  # because the automatic detection only recognises `depot` / `central`.
  dosing <- c("depot", "depot_osimertinib")

  # Paper-mechanistic QSP states with no canonical role in
  # inst/references/compartment-names.md. `egfr_signal` (rather than `egfr`)
  # because the canonical `egfr` is the renal estimated-glomerular-filtration-
  # rate state; this state is the tumour EGFR signalling level relative to its
  # untreated baseline. The four `*_unperturbed` states are the drug-free twin
  # tumour chain of Hu 2026 supplementary Equations S1-S4, carried solely so
  # that real-time TGI% (Equation S6) is available as a model output.
  paper_specific_compartments <- c(
    "enzyme_folate", "folate", "egfr_signal", "total_death",
    "cycling_cells_unperturbed", "damaged_cells1_unperturbed",
    "damaged_cells2_unperturbed", "damaged_cells3_unperturbed"
  )

  units <- list(
    time          = "day",
    dosing        = "mg/kg (pemetrexed into depot; osimertinib into depot_osimertinib)",
    concentration = "mg/L (pemetrexed, Cc); ug/L (osimertinib, Cc_osimertinib); the tumour output tumor_vol is a volume in mm3"
  )

  compartmentData <- list(
    depot                      = list(analyte = "pemetrexed",  units = "mg/kg", specimen = "administration site", verified = TRUE),
    central                    = list(analyte = "pemetrexed",  units = "mg/kg", specimen = "plasma",              verified = TRUE),
    peripheral1                = list(analyte = "pemetrexed",  units = "mg/kg", specimen = "plasma",              verified = TRUE),
    depot_osimertinib          = list(analyte = "osimertinib", units = "mg/kg", specimen = "administration site", verified = TRUE),
    central_osimertinib        = list(analyte = "osimertinib", units = "mg/kg", specimen = "plasma",              verified = TRUE),
    peripheral1_osimertinib    = list(analyte = "osimertinib", units = "mg/kg", specimen = "plasma",              verified = TRUE),
    enzyme_folate              = list(analyte = "folate-metabolizing enzymes", units = "unitless (fraction of untreated baseline)", specimen = "tumor", verified = TRUE),
    folate                     = list(analyte = "folate",      units = "unitless (fraction of untreated baseline)", specimen = "tumor", verified = TRUE),
    egfr_signal                = list(analyte = "EGFR",        units = "unitless (fraction of untreated baseline)", specimen = "tumor", verified = TRUE),
    cycling_cells              = list(analyte = "cells",       units = "mm3",   specimen = "tumor",                verified = TRUE),
    damaged_cells1             = list(analyte = "cells",       units = "mm3",   specimen = "tumor",                verified = TRUE),
    damaged_cells2             = list(analyte = "cells",       units = "mm3",   specimen = "tumor",                verified = TRUE),
    damaged_cells3             = list(analyte = "cells",       units = "mm3",   specimen = "tumor",                verified = TRUE),
    total_death                = list(analyte = "cells",       units = "mm3",   specimen = "not applicable",       verified = TRUE),
    cycling_cells_unperturbed  = list(analyte = "cells",       units = "mm3",   specimen = "not applicable",       verified = TRUE),
    damaged_cells1_unperturbed = list(analyte = "cells",       units = "mm3",   specimen = "not applicable",       verified = TRUE),
    damaged_cells2_unperturbed = list(analyte = "cells",       units = "mm3",   specimen = "not applicable",       verified = TRUE),
    damaged_cells3_unperturbed = list(analyte = "cells",       units = "mm3",   specimen = "not applicable",       verified = TRUE)
  )

  # Hu 2026 builds one model with one parameter set that reproduces every
  # treatment arm; no covariate enters the structural model. The BIM-deletion
  # polymorphism and the reduced-OSI-sensitivity scenarios of Figures 9 and 11
  # are parameter-perturbation SIMULATION scenarios (typical values of kbim /
  # gamma_bim multiplied by 0.1 / 5, EC50_osi by 5, Imax_osi by 1/5), not
  # covariate effects fitted in the model -- see the vignette.
  covariateData <- list()

  population <- list(
    species        = "mouse (female BALB/c nude with subcutaneous HCC827 EGFR-mutant NSCLC xenograft); PK parameters from PC9-bearing male/female BALB/c nude mice",
    n_subjects     = 25L,
    n_studies      = 1L,
    age_range      = "7 weeks at inoculation",
    weight_range   = "approximately 20 g",
    sex_female_pct = 100,
    race_ethnicity = NA,
    disease_state  = "subcutaneous HCC827 human EGFR-exon-19-deletion NSCLC xenograft (1e7 cells in 50% high-concentration Matrigel, right flank); randomised at a tumour volume of approximately 200 mm3",
    dose_range     = "pemetrexed 35 mg/kg intraperitoneally three times daily 4 h apart (105 mg/kg/day) on days 0-1 of each 7-day cycle; osimertinib 1 mg/kg orally once daily on days 0-2 (concurrent) or days 2-4 (48 h sequential) of each cycle; three cycles",
    regions        = "preclinical (China Pharmaceutical University, Nanjing, China; animal ethics approval YSL-202504062)",
    notes          = "Five arms of n = 5 (control, PEM, OSI, PEM + OSI, PEM -> OSI); tumour volume V = (pi/6) * a * b^2 measured by caliper every 3 days to day 18 (Hu 2026 Figure 6B, replotted from the authors' earlier report doi:10.1016/j.canlet.2024.217124). The PK sub-models were fitted separately in WinNonlin to pooled mean plasma profiles from PC9-bearing BALB/c nude mice dosed pemetrexed 100 mg/kg i.p. and osimertinib 5 mg/kg p.o. (Hu 2026 Figures 5D and S4A,B). Between-subject variability is the 30% CV log-normal spread applied to every PD and tumour-growth parameter in the Monte Carlo analysis of Hu 2026 Supplementary Method S1.7."
  )

  ini({
    # ------------------------------------------------------------------
    # Fixed-versus-estimated follows the "Source" column of Hu 2026
    # Table S19: rows sourced from "Assumed Value" or "Literature" are
    # wrapped in fixed(); rows sourced from "Curve Fitting", "PK Analysis",
    # "in vitro data" or "IVIVE" are left free.
    #
    # Every PD / tumour-growth parameter is log-transformed so that the
    # 30% CV log-normal between-subject variability of Hu 2026
    # Supplementary Method S1.7 is a mu-referenced eta.
    # ------------------------------------------------------------------

    # --- Pemetrexed PK (two-compartment, intraperitoneal, F fixed at 1) ----
    # Hu 2026 Equations (1)-(3); values from Table S19 "PEM PK Module".
    lka  <- log(102.3379); label("Pemetrexed intraperitoneal absorption rate ka (1/day; log-transformed)")            # Table S19 ka_pem = 102.3379 1/day (PK Analysis)
    lkel <- log(56.6998);  label("Pemetrexed elimination rate from central kel (1/day; log-transformed)")             # Table S19 kel_pem = 56.6998 1/day (PK Analysis)
    lk12 <- log(2.3398);   label("Pemetrexed central-to-peripheral rate k12 (1/day; log-transformed)")                # Table S19 k12_pem = 2.3398 1/day (PK Analysis)
    lk21 <- log(2.8223);   label("Pemetrexed peripheral-to-central rate k21 (1/day; log-transformed)")                # Table S19 k21_pem = 2.8223 1/day (PK Analysis)
    lvc  <- log(0.39125);  label("Pemetrexed central volume of distribution (L/kg; log-transformed)")                 # Table S19 PEM_V1 = 0.39125 L/kg (PK Analysis)
    lvp  <- log(0.32436);  label("Pemetrexed peripheral volume of distribution (L/kg; log-transformed)")              # Table S19 PEM_V2 = 0.32436 L/kg (PK Analysis)

    # --- Osimertinib PK (two-compartment, oral; volumes are apparent V/F) ---
    # Hu 2026 Equations (4)-(6); values from Table S19 "OSI PK Module".
    lka_osimertinib  <- log(24.6514); label("Osimertinib oral absorption rate ka (1/day; log-transformed)")           # Table S19 ka_osi = 24.6514 1/day (PK Analysis)
    lkel_osimertinib <- log(10.6071); label("Osimertinib elimination rate from central kel (1/day; log-transformed)") # Table S19 k_el_osi = 10.6071 1/day (PK Analysis)
    lk12_osimertinib <- log(11.005);  label("Osimertinib central-to-peripheral rate k12 (1/day; log-transformed)")    # Table S19 k12_osi = 11.005 1/day (PK Analysis)
    lk21_osimertinib <- log(5.4194);  label("Osimertinib peripheral-to-central rate k21 (1/day; log-transformed)")    # Table S19 k21_osi = 5.4194 1/day (PK Analysis)
    lvc_osimertinib  <- log(7.1758);  label("Osimertinib apparent central volume V1/F (L/kg; log-transformed)")       # Table S19 OSI_V1/F = 7.1758 L/kg (PK Analysis)
    lvp_osimertinib  <- log(14.5716); label("Osimertinib apparent peripheral volume V2/F (L/kg; log-transformed)")    # Table S19 OSI_V2/F = 14.5716 L/kg (PK Analysis)

    # --- Folate module (Hu 2026 Equations (7)-(8); Table S19 "Folate Module") ---
    lkout_enzyme  <- fixed(log(1));       label("Turnover rate of the folate-metabolizing enzyme pool (1/day; log-transformed)")        # Table S19 kout_Enzyme = 1 1/day (Assumed Value)
    lkout_folate  <- fixed(log(3));       label("Turnover rate of the folate pool (1/day; log-transformed)")                           # Table S19 kout_folate = 3 1/day (Assumed Value)
    lgamma_enzyme <- log(2.0129);         label("Power of enzyme level on folate synthesis (unitless; log-transformed)")               # Table S19 gamma_Enzyme = 2.0129 (IVIVE)
    lemax_enzyme  <- log(5.56);           label("Pemetrexed Emax on folate-enzyme degradation (unitless; log-transformed)")            # Table S19 Emax_pem = 5.56 (IVIVE)
    lec50_enzyme  <- log(0.47315);        label("Pemetrexed EC50 for folate-enzyme degradation (mg/L plasma; log-transformed)")        # Table S19 EC50_pem,plasma = 0.47315 mg/L (IVIVE from EC50_pem,medium = 300 nM; Supplementary Method S1.6)

    # --- EGFR module (Hu 2026 Equation (9); Table S19 "EGFR Module") ---
    lkout_egfr_signal           <- log(1.5);          label("Turnover rate of the tumour EGFR signal (1/day; log-transformed)")                          # Table S19 kout_EGFR = 1.5 1/day (IVIVE, Figure 7E,F)
    lkfeedback_egfr_signal      <- log(0.2);          label("Damaged-cell-driven fractional increase in EGFR synthesis (unitless; log-transformed)")     # Table S19 k_EGFR_feedback = 0.2 (Curve Fitting)
    lgamma_feedback_egfr_signal <- log(1.45);         label("Power on the damaged-cell fraction in the EGFR feedback term (unitless; log-transformed)")  # Table S19 gamma_EGFR_feedback = 1.45 (Curve Fitting)
    limax_osimertinib           <- fixed(log(22.872)); label("Osimertinib Imax on EGFR degradation (unitless; log-transformed)")                         # Table S19 Imax_osi = 22.872 (Literature, doi:10.1158/1535-7163.MCT-16-0142)
    lec50_osimertinib           <- log(48.866);       label("Osimertinib EC50 for EGFR inhibition (ug/L plasma; log-transformed)")                       # Table S19 EC50_osi,plasma = 48.866 ug/L (IVIVE from EC50_osi,medium = 14.49 nM)
    lgamma_osimertinib          <- fixed(log(2));     label("Hill coefficient on the osimertinib EGFR-inhibition term (unitless; log-transformed)")      # Table S19 gamma_osi = 2 (Assumed Value)

    # --- Tumour-growth-inhibition module (Hu 2026 Equations (10)-(17); Table S19 "TGI Module") ---
    ltumorExpGrowth           <- log(0.1032);  label("Simeoni exponential-phase tumour growth rate lambda0 (1/day; log-transformed)")            # Table S19 Lambda_0 = 0.1032 1/day (Curve Fitting, Figure 8E)
    ltumorLinGrowth           <- log(51.08);   label("Simeoni linear-phase tumour growth rate lambda1 (mm3/day; log-transformed)")               # Table S19 Lambda_1 = 51.08 mm3/day (Curve Fitting, Figure 8E)
    ldamageRate               <- log(0.016);   label("Natural proliferating-to-damaged transition rate k1 (1/day; log-transformed)")             # Table S19 k1 = 0.016 1/day (Curve Fitting)
    ldamageTransit            <- log(0.0045);  label("Natural damaged-cell transit / apoptosis rate k2 (1/day; log-transformed)")                # Table S19 k2 = 0.0045 1/day (Curve Fitting)
    lemax_folate              <- log(49.5);    label("Maximum fractional increase in k1 from folate depletion (unitless; log-transformed)")      # Table S19 Emax_folate = 49.5 (IVIVE)
    lec50_folate              <- fixed(log(0.5));  label("EC50 of the folate-depletion effect on k1 (fraction below baseline; log-transformed)") # Table S19 EC50_folate = 0.5 (Assumed Value)
    lkbim                     <- fixed(log(211));  label("Maximum fractional increase in k2 from EGFR depletion (Bim surrogate; log-transformed)") # Table S19 k_bim = 211 (Literature, doi:10.1158/1535-7163.MCT-16-0142 f_kill)
    lgamma_bim                <- log(0.595);   label("Power on (1 - EGFR) in the Bim-driven apoptosis term (unitless; log-transformed)")         # Table S19 gamma_bim = 0.595 (IVIVE, Figure 7L)
    lgamma_g1_egfr_signal     <- fixed(log(8)); label("Power of EGFR level on the k1 transition (G1-arrest antagonism; log-transformed)")        # Table S19 gamma_G1 = 8 (Assumed Value)
    lgamma_prolif_egfr_signal <- log(3.998);   label("Power of EGFR level on the proliferation rate of X1 (unitless; log-transformed)")          # Table S19 gamma_EGFR = 3.998 (IVIVE, Figure 7A,B)
    psi                       <- fixed(20);    label("Simeoni switching exponent psi (unitless)")                                                # Table S19 psi = 20 (Literature, Simeoni 2004 doi:10.1158/0008-5472.CAN-03-2524)
    tumor_vol0                <- fixed(200);   label("Tumour volume at randomisation (mm3)")                                                     # Table S20 initial value X1 = 200 mm3

    # --- Between-subject variability -------------------------------------
    # Hu 2026 Supplementary Method S1.7: every PD and tumour-growth-related
    # parameter was given 30% CV log-normal inter-individual variability,
    # with sigma^2 = ln(1 + CV^2) = ln(1.09) = 0.0861777 (Equations S10-S12).
    # The variance is fixed() because the paper prescribes 30% CV rather than
    # estimating it. The switching exponent psi and the initial tumour volume
    # carry no eta: psi is a numerical smoothing exponent taken from Simeoni
    # 2004 rather than a biological PD parameter, and Method S1.7 does not
    # list the baseline tumour volume among the varied quantities.
    # No eta is placed on the PK parameters -- the mouse PK was fitted to
    # pooled mean profiles, so no individual PK variability is reported.
    etalgamma_enzyme               ~ fixed(0.0861777)  # Supplementary Method S1.7: 30% CV log-normal
    etalemax_enzyme                ~ fixed(0.0861777)  # Supplementary Method S1.7: 30% CV log-normal
    etalec50_enzyme                ~ fixed(0.0861777)  # Supplementary Method S1.7: 30% CV log-normal
    etalkout_enzyme                ~ fixed(0.0861777)  # Supplementary Method S1.7: 30% CV log-normal
    etalkout_folate                ~ fixed(0.0861777)  # Supplementary Method S1.7: 30% CV log-normal
    etalkout_egfr_signal           ~ fixed(0.0861777)  # Supplementary Method S1.7: 30% CV log-normal
    etalkfeedback_egfr_signal      ~ fixed(0.0861777)  # Supplementary Method S1.7: 30% CV log-normal
    etalgamma_feedback_egfr_signal ~ fixed(0.0861777)  # Supplementary Method S1.7: 30% CV log-normal
    etalimax_osimertinib           ~ fixed(0.0861777)  # Supplementary Method S1.7: 30% CV log-normal
    etalec50_osimertinib           ~ fixed(0.0861777)  # Supplementary Method S1.7: 30% CV log-normal
    etalgamma_osimertinib          ~ fixed(0.0861777)  # Supplementary Method S1.7: 30% CV log-normal
    etaltumorExpGrowth             ~ fixed(0.0861777)  # Supplementary Method S1.7: 30% CV log-normal
    etaltumorLinGrowth             ~ fixed(0.0861777)  # Supplementary Method S1.7: 30% CV log-normal
    etaldamageRate                 ~ fixed(0.0861777)  # Supplementary Method S1.7: 30% CV log-normal
    etaldamageTransit              ~ fixed(0.0861777)  # Supplementary Method S1.7: 30% CV log-normal
    etalemax_folate                ~ fixed(0.0861777)  # Supplementary Method S1.7: 30% CV log-normal
    etalec50_folate                ~ fixed(0.0861777)  # Supplementary Method S1.7: 30% CV log-normal
    etalkbim                       ~ fixed(0.0861777)  # Supplementary Method S1.7: 30% CV log-normal (two-point BIM-deletion mixture applied on top; see vignette)
    etalgamma_bim                  ~ fixed(0.0861777)  # Supplementary Method S1.7: 30% CV log-normal (two-point BIM-deletion mixture applied on top; see vignette)
    etalgamma_g1_egfr_signal       ~ fixed(0.0861777)  # Supplementary Method S1.7: 30% CV log-normal
    etalgamma_prolif_egfr_signal   ~ fixed(0.0861777)  # Supplementary Method S1.7: 30% CV log-normal

    # --- Residual error ---------------------------------------------------
    propSd_tumor_vol <- fixed(0); label("Proportional residual error on tumour volume (fraction; ZERO - not reported in source)")  # Hu 2026 reports no residual-error model: the QSP model is calibrated to mean tumour volumes and used for simulation only
  })

  model({
    # ==================================================================
    # 1. Individual parameters
    # ==================================================================
    ka  <- exp(lka)
    kel <- exp(lkel)
    k12 <- exp(lk12)
    k21 <- exp(lk21)
    vc  <- exp(lvc)
    vp  <- exp(lvp)

    ka_osimertinib  <- exp(lka_osimertinib)
    kel_osimertinib <- exp(lkel_osimertinib)
    k12_osimertinib <- exp(lk12_osimertinib)
    k21_osimertinib <- exp(lk21_osimertinib)
    vc_osimertinib  <- exp(lvc_osimertinib)
    vp_osimertinib  <- exp(lvp_osimertinib)

    kout_enzyme  <- exp(lkout_enzyme  + etalkout_enzyme)
    kout_folate  <- exp(lkout_folate  + etalkout_folate)
    gamma_enzyme <- exp(lgamma_enzyme + etalgamma_enzyme)
    emax_enzyme  <- exp(lemax_enzyme  + etalemax_enzyme)
    ec50_enzyme  <- exp(lec50_enzyme  + etalec50_enzyme)

    kout_egfr_signal           <- exp(lkout_egfr_signal           + etalkout_egfr_signal)
    kfeedback_egfr_signal      <- exp(lkfeedback_egfr_signal      + etalkfeedback_egfr_signal)
    gamma_feedback_egfr_signal <- exp(lgamma_feedback_egfr_signal + etalgamma_feedback_egfr_signal)
    imax_osimertinib           <- exp(limax_osimertinib           + etalimax_osimertinib)
    ec50_osimertinib           <- exp(lec50_osimertinib           + etalec50_osimertinib)
    gamma_osimertinib          <- exp(lgamma_osimertinib          + etalgamma_osimertinib)

    tumorExpGrowth           <- exp(ltumorExpGrowth           + etaltumorExpGrowth)
    tumorLinGrowth           <- exp(ltumorLinGrowth           + etaltumorLinGrowth)
    damageRate               <- exp(ldamageRate               + etaldamageRate)
    damageTransit            <- exp(ldamageTransit            + etaldamageTransit)
    emax_folate              <- exp(lemax_folate              + etalemax_folate)
    ec50_folate              <- exp(lec50_folate              + etalec50_folate)
    kbim                     <- exp(lkbim                     + etalkbim)
    gamma_bim                <- exp(lgamma_bim                + etalgamma_bim)
    gamma_g1_egfr_signal     <- exp(lgamma_g1_egfr_signal     + etalgamma_g1_egfr_signal)
    gamma_prolif_egfr_signal <- exp(lgamma_prolif_egfr_signal + etalgamma_prolif_egfr_signal)

    # ==================================================================
    # 2. Pemetrexed PK -- Hu 2026 Equations (1)-(3)
    #    States are amounts per kg (mg/kg); bioavailability is fixed at 1
    #    for the intraperitoneal route (Hu 2026 Section 2.3.2).
    # ==================================================================
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - k12 * central + k21 * peripheral1 - kel * central
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    Cc <- central / vc                        # mg/L, Hu 2026 Section 2.3.2: Cpem,1 = Xpem,1 / V1,pem

    # ==================================================================
    # 3. Osimertinib PK -- Hu 2026 Equations (4)-(6)
    #    Oral bioavailability is absorbed into the apparent volumes V/F
    #    (Table S19 reports OSI_V1/F and OSI_V2/F), matching the ODE export
    #    in Table S17 which carries no explicit Fa term.
    # ==================================================================
    d/dt(depot_osimertinib)       <- -ka_osimertinib * depot_osimertinib
    d/dt(central_osimertinib)     <-  ka_osimertinib * depot_osimertinib -
                                      k12_osimertinib * central_osimertinib +
                                      k21_osimertinib * peripheral1_osimertinib -
                                      kel_osimertinib * central_osimertinib
    d/dt(peripheral1_osimertinib) <-  k12_osimertinib * central_osimertinib -
                                      k21_osimertinib * peripheral1_osimertinib
    # 1000 converts mg/L to ug/L so that Cc_osimertinib is on the same scale
    # as EC50_osi,plasma = 48.866 ug/L (Table S19).
    Cc_osimertinib <- 1000 * central_osimertinib / vc_osimertinib   # ug/L

    # ==================================================================
    # 4. Folate module -- Hu 2026 Equations (7)-(8), Table S17 rows 13-14
    #    Both states are levels relative to the untreated baseline of 1.
    #    gamma_pem is fixed at 1 and therefore absent from the Hill term
    #    (footnote to Tables S17 and S21).
    # ==================================================================
    d/dt(enzyme_folate) <- kout_enzyme -
      kout_enzyme * enzyme_folate * (1 + emax_enzyme * Cc / (ec50_enzyme + Cc))
    d/dt(folate) <- kout_folate * enzyme_folate^gamma_enzyme - kout_folate * folate

    # ==================================================================
    # 5. Tumour bookkeeping -- Hu 2026 Equations (14)-(16), Table S18
    #    Damaged_percent is the FRACTION X_damaged / X_total (Table S18 row 3
    #    and Table S20, which give it as unitless with initial value 0); the
    #    "x 100%" in printed Equation (15) is cosmetic. Feeding a percentage
    #    into the feedback term with k_EGFR_feedback = 0.2 and
    #    gamma_EGFR_feedback = 1.45 would raise baseline EGFR a hundredfold.
    # ==================================================================
    #
    # Numerical guard on the denominator: in long simulations (the 210-day
    # Monte Carlo of Hu 2026 Section 3.6) the most responsive virtual animals
    # drive the tumour to ~1e-24 mm3, where solver round-off pushes the states
    # through zero and the published ratio becomes 0/0 (and the subsequent
    # damaged_frac^gamma_feedback becomes NaN for a negative base). A single
    # tumour cell is of order 1e-6 mm3, so flooring the denominator at
    # 1e-12 mm3 and the numerator at zero leaves the ratio unchanged at every
    # physically meaningful tumour volume.
    x_total       <- cycling_cells + damaged_cells1 + damaged_cells2 + damaged_cells3
    x_damaged     <- damaged_cells1 + damaged_cells2 + damaged_cells3
    x_total_pos   <- max(x_total, 1e-12)
    damaged_frac  <- max(x_damaged, 0) / x_total_pos

    # ==================================================================
    # 6. EGFR module -- Hu 2026 Equation (9), Table S17 row 18
    # ==================================================================
    d/dt(egfr_signal) <-
      kout_egfr_signal * (1 + kfeedback_egfr_signal * damaged_frac^gamma_feedback_egfr_signal) -
      kout_egfr_signal * egfr_signal *
        (1 + imax_osimertinib * Cc_osimertinib^gamma_osimertinib /
               (ec50_osimertinib^gamma_osimertinib + Cc_osimertinib^gamma_osimertinib))

    # ==================================================================
    # 7. Tumour growth inhibition -- Hu 2026 Equations (10)-(13) and (17),
    #    Table S17 rows 1-4 and 9.
    # ==================================================================
    # Simeoni switching growth term. Hu 2026 writes the switch on X1 (the
    # proliferating pool) rather than on the total tumour volume; Table S17
    # row 1 and printed Equation (10) agree on this, so X1 is used here.
    prolif <- tumorExpGrowth * cycling_cells /
      (1 + (tumorExpGrowth * cycling_cells / tumorLinGrowth)^psi)^(1 / psi)

    # Folate-depletion stimulation of the X1 -> X2 transition.
    # gamma_folate is fixed at 1 and therefore absent (footnote to Table S17).
    folate_effect <- 1 + emax_folate * (1 - folate) / (ec50_folate + (1 - folate))

    # Bim-driven acceleration of transit through the damaged-cell chain.
    # Hu 2026 uses 1 - EGFR as a surrogate for Bim (Section 2.3.6, Figure 7J);
    # the power law gamma_bim = 0.595 was fitted on Figure 7L over EGFR levels
    # BELOW baseline. Damaged-cell feedback drives EGFR above 1 whenever
    # osimertinib is absent, where (1 - EGFR)^0.595 has no real value, so the
    # deviation is clamped at zero: above baseline there is no Bim induction
    # and hence no acceleration. Without the clamp the published equation
    # returns NaN within the first solver step of every arm, including control.
    egfr_deviation <- max(1 - egfr_signal, 0)
    bim_effect     <- 1 + kbim * egfr_deviation^gamma_bim

    # X1 -> X2 flux: folate depletion accelerates it, EGFR-driven G1 arrest
    # (egfr_signal^gamma_g1) antagonises it.
    damage_flux <- damageRate * cycling_cells * folate_effect *
      egfr_signal^gamma_g1_egfr_signal

    d/dt(cycling_cells)  <- prolif * egfr_signal^gamma_prolif_egfr_signal - damage_flux
    d/dt(damaged_cells1) <- damage_flux - damageTransit * damaged_cells1 * bim_effect
    d/dt(damaged_cells2) <- damageTransit * bim_effect * (damaged_cells1 - damaged_cells2)
    d/dt(damaged_cells3) <- damageTransit * bim_effect * (damaged_cells2 - damaged_cells3)
    d/dt(total_death)    <- damageTransit * damaged_cells3 * bim_effect

    # ==================================================================
    # 8. Unperturbed (drug-free) twin chain -- Hu 2026 Supplementary
    #    Equations (S1)-(S6), Table S17 rows 5-8. Identical parameters,
    #    no drug effects; used only to compute real-time TGI%.
    # ==================================================================
    prolif_unperturbed <- tumorExpGrowth * cycling_cells_unperturbed /
      (1 + (tumorExpGrowth * cycling_cells_unperturbed / tumorLinGrowth)^psi)^(1 / psi)

    d/dt(cycling_cells_unperturbed)  <- prolif_unperturbed - damageRate * cycling_cells_unperturbed
    d/dt(damaged_cells1_unperturbed) <- damageRate * cycling_cells_unperturbed -
                                        damageTransit * damaged_cells1_unperturbed
    d/dt(damaged_cells2_unperturbed) <- damageTransit *
                                        (damaged_cells1_unperturbed - damaged_cells2_unperturbed)
    d/dt(damaged_cells3_unperturbed) <- damageTransit *
                                        (damaged_cells2_unperturbed - damaged_cells3_unperturbed)

    x_total_unperturbed <- cycling_cells_unperturbed + damaged_cells1_unperturbed +
                           damaged_cells2_unperturbed + damaged_cells3_unperturbed

    # ==================================================================
    # 9. Initial conditions -- Hu 2026 Table S20
    # ==================================================================
    cycling_cells(0)              <- tumor_vol0
    cycling_cells_unperturbed(0)  <- tumor_vol0
    enzyme_folate(0)              <- 1
    folate(0)                     <- 1
    egfr_signal(0)                <- 1

    # ==================================================================
    # 10. Derived outputs and observation
    # ==================================================================
    tgi         <- 1 - x_total / x_total_unperturbed   # Table S18 row 5 (fraction)
    x1_fraction <- max(cycling_cells, 0) / x_total_pos # Table S18 row 6 (fraction; same guard)
    tumor_vol   <- x_total                             # mm3, the observed endpoint (Figure 6B)

    tumor_vol ~ prop(propSd_tumor_vol)
  })
}
