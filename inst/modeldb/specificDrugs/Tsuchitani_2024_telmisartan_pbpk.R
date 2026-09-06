Tsuchitani_2024_telmisartan_pbpk <- function() {
  description <- paste(
    "PBPK-TMDD (semi-mechanistic, five-compartment tandem liver,",
    "segregated-flow intestine). Telmisartan and its 1-O-acylglucuronide",
    "(Tel-GLU) after single oral (5-160 mg) and intravenous (10-120 mg)",
    "doses in healthy volunteers, fitted top-down and then middle-out with",
    "the Cluster Gauss-Newton method (CGNM). The model explains the",
    "nonlinear pharmacokinetics of telmisartan by two saturable processes.",
    "First, hepatic uptake: sinusoid-to-hepatocyte transport is",
    "Michaelis-Menten OATP1B3 uptake plus passive diffusion, and the in",
    "vivo Km,OATP1B3 is estimated at 3.51 nM, far below the 810 nM measured",
    "in vitro in 0.3 percent albumin, which the authors attribute to the",
    "albumin-mediated uptake effect and confirm with a plated-human-",
    "hepatocyte Km of 4.49 nM in 4.5 percent albumin. Second, target-",
    "mediated drug disposition: a finite pool of AT1 receptors (R_total",
    "13.4 umol) is distributed over the central compartment (57.2 percent)",
    "and the blood-accessible tissues, and reversible binding with Kd",
    "0.160 nM saturates over the therapeutic range, driving the",
    "nonlinearity below 20 mg. The liver is five sinusoidal",
    "sub-compartments in series, each exchanging with its own hepatocyte",
    "sub-compartment; the intestine is a segregated-flow model with",
    "duodenum, jejunum and ileum lumen, enterocyte and mucosal-blood",
    "compartments plus caecum and colon; muscle, skin and adipose are",
    "perfusion-limited. Enterohepatic circulation of both telmisartan and",
    "Tel-GLU is included, with deconjugation of Tel-GLU in the hepatocytes",
    "and by gut microbiota. The 17 parameters without fixed() are the",
    "middle-out CGNM medians of main-text Table 1; every other value is a",
    "fixed physiological or compound constant from Tables S3-S7. CGNM is a",
    "fixed-effects nonlinear-least-squares method, so this model has no",
    "between-subject variability and the paper reports no residual-error",
    "model -- the propSd term is a placeholder. See the vignette Errata."
  )
  reference <- paste(
    "Tsuchitani T, Tomaru A, Aoki Y, Ishiguro N, Tsuda Y, Sugiyama Y.",
    "Elucidating nonlinear pharmacokinetics of telmisartan: Integration of",
    "target-mediated drug disposition and OATP1B3-mediated hepatic uptake",
    "in a physiologically based model.",
    "CPT Pharmacometrics Syst Pharmacol. 2024;13(7):1224-1237.",
    "doi:10.1002/psp4.13154.",
    "The ODE system and the auxiliary parameter definitions are",
    "transcribed from Data S1 (Supporting Information file",
    "PSP4-13-1224-s002.docx), which carries the authors' complete model",
    "code in d/dt() form. The estimated parameters are the 'Final",
    "parameters median' column of main-text Table 1; the fixed",
    "physiological and compound constants are Supplementary Tables S3-S7",
    "in Appendix S1 (PSP4-13-1224-s001.docx).",
    sep = " "
  )
  vignette <- "Tsuchitani_2024_telmisartan_pbpk"

  # Data S1 numbers its states y1..y35 (plus a "G" infix for the
  # Tel-GLU species and R_complex / R_free / RO_ prefixes for the
  # receptor layer). The mapping to the names used below is given in
  # the comment on every d/dt() line. The stems follow the canonical
  # PBPK sub-compartment vocabulary already used by this laboratory's
  # two registered models (Aoki_2024_bosentan_pbpk.R,
  # Tsuchitani_2026_apixaban_pbpk.R): `is_liver<n>` is the hepatic
  # extracellular (sinusoidal) space and `int_liver<n>` the
  # intracellular (hepatocyte) space of dispersion segment n; only the
  # 1..5 index is paper-specific. `_gluc` is the registered
  # glucuronide-metabolite suffix (R/conventions.R
  # registeredMetabolites) and carries Tel-GLU.
  #
  # The receptor layer uses the `target_<location>` / `complex_<location>`
  # TMDD convention, but over eleven anatomical locations rather than the
  # two (`csf`, `isf`) that `targetLocationRegex` in R/conventions.R
  # enumerates, so the states are declared paper-specific here rather
  # than widening that register for a single paper. `occupancy_<location>`
  # is fractional receptor occupancy, carried as an explicit state
  # exactly as Data S1 carries it (`RO_*`); it is mathematically
  # redundant with complex_<location> / (r_total * at1_<location>) and
  # the vignette asserts that identity as a numerical check on the
  # integration.
  paper_specific_compartments <- c(
    "is_liver1", "is_liver2", "is_liver3", "is_liver4", "is_liver5",
    "int_liver1", "int_liver2", "int_liver3", "int_liver4", "int_liver5",
    "int_liver1_gluc", "int_liver2_gluc", "int_liver3_gluc",
    "int_liver4_gluc", "int_liver5_gluc",
    "duodenum_lumen", "jejunum_lumen", "ileum_lumen",
    "caecum_lumen", "colon_lumen",
    "duodenum_lumen_gluc", "jejunum_lumen_gluc", "ileum_lumen_gluc",
    "caecum_lumen_gluc", "colon_lumen_gluc",
    "duodenum_ent", "jejunum_ent", "ileum_ent",
    "duodenum_ent_gluc", "jejunum_ent_gluc", "ileum_ent_gluc",
    "duodenum_muc", "jejunum_muc", "ileum_muc",
    "serosa",
    "ehc1", "ehc2", "ehc3",
    "ehc1_gluc", "ehc2_gluc", "ehc3_gluc",
    "a_feces", "a_feces_gluc", "a_urine", "auc_blood",
    "target_central", "target_liver1", "target_liver2", "target_liver3",
    "target_liver4", "target_liver5", "target_duodenum",
    "target_jejunum", "target_ileum", "target_muscle", "target_skin",
    "target_adipose",
    "complex_central", "complex_liver1", "complex_liver2",
    "complex_liver3", "complex_liver4", "complex_liver5",
    "complex_duodenum", "complex_jejunum", "complex_ileum",
    "complex_muscle", "complex_skin", "complex_adipose",
    "occupancy_central", "occupancy_liver1", "occupancy_liver2",
    "occupancy_liver3", "occupancy_liver4", "occupancy_liver5",
    "occupancy_duodenum", "occupancy_jejunum", "occupancy_ileum",
    "occupancy_muscle", "occupancy_skin", "occupancy_adipose"
  )

  # Time in hours. Data S1 works in molar units throughout: the tissue,
  # lumen and enterocyte states hold concentrations in umol/L
  # (micromolar) and the receptor, bile, urine and faecal states hold
  # amounts in umol, so doses must be supplied in umol. Telmisartan
  # MW 514.6 g/mol, so mg * 1000 / 514.6 = umol; a 40 mg tablet is
  # 77.7 umol. That molecular weight is confirmed twice by the paper
  # itself: the LC-MS/MS precursor ion is 515.143 for [M+H]+
  # (Materials and Methods, "Quantification of telmisartan by
  # LC-MS/MS"), and Table S1 reports CLh = Dose/AUCinf = 49.5 L/h for
  # the 10 mg i.v. arm against AUCinf = 393 nM*h, which requires a dose
  # of 19.45 umol, i.e. 10 mg / 514.1 g/mol.
  units <- list(time = "h", dosing = "umol", concentration = "umol/L")

  # No covariates. Every physiological volume, flow and partition
  # coefficient in Tables S3-S7 is a fixed constant for a single typical
  # 78 kg adult; the published model carries no body-weight or
  # demographic scaling of any kind. The pharmacogenomic simulation of
  # Figure 6 is performed by scaling vmax_oatp1b3 / vmax_ugt outside the
  # model rather than by a genotype covariate in it (see the vignette).
  covariateData <- list()

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. analyte/specimen proposed by a local model from
  # the model description; units derived from the units block.
  # verified = FALSE means NOT checked against the source paper.
  compartmentData <- list(
    central             = list(analyte = "telmisartan", units = "umol/L", specimen = "whole blood", verified = FALSE),
    is_liver1           = list(analyte = "telmisartan", units = "umol/L", specimen = "tissue", verified = FALSE),
    is_liver2           = list(analyte = "telmisartan", units = "umol/L", specimen = "tissue", verified = FALSE),
    is_liver3           = list(analyte = "telmisartan", units = "umol/L", specimen = "tissue", verified = FALSE),
    is_liver4           = list(analyte = "telmisartan", units = "umol/L", specimen = "tissue", verified = FALSE),
    is_liver5           = list(analyte = "telmisartan", units = "umol/L", specimen = "tissue", verified = FALSE),
    int_liver1          = list(analyte = "telmisartan", units = "umol/L", specimen = "tissue", verified = FALSE),
    int_liver2          = list(analyte = "telmisartan", units = "umol/L", specimen = "tissue", verified = FALSE),
    int_liver3          = list(analyte = "telmisartan", units = "umol/L", specimen = "tissue", verified = FALSE),
    int_liver4          = list(analyte = "telmisartan", units = "umol/L", specimen = "tissue", verified = FALSE),
    int_liver5          = list(analyte = "telmisartan", units = "umol/L", specimen = "tissue", verified = FALSE),
    int_liver1_gluc     = list(analyte = "telmisartan glucuronide", units = "umol/L", specimen = "tissue", verified = FALSE),
    int_liver2_gluc     = list(analyte = "telmisartan glucuronide", units = "umol/L", specimen = "tissue", verified = FALSE),
    int_liver3_gluc     = list(analyte = "telmisartan glucuronide", units = "umol/L", specimen = "tissue", verified = FALSE),
    int_liver4_gluc     = list(analyte = "telmisartan glucuronide", units = "umol/L", specimen = "tissue", verified = FALSE),
    int_liver5_gluc     = list(analyte = "telmisartan glucuronide", units = "umol/L", specimen = "tissue", verified = FALSE),
    duodenum_lumen      = list(analyte = "telmisartan", units = "umol/L", specimen = "administration site", verified = FALSE),
    jejunum_lumen       = list(analyte = "telmisartan", units = "umol/L", specimen = "administration site", verified = FALSE),
    ileum_lumen         = list(analyte = "telmisartan", units = "umol/L", specimen = "administration site", verified = FALSE),
    caecum_lumen        = list(analyte = "telmisartan", units = "umol/L", specimen = "administration site", verified = FALSE),
    colon_lumen         = list(analyte = "telmisartan", units = "umol/L", specimen = "administration site", verified = FALSE),
    duodenum_lumen_gluc = list(analyte = "telmisartan glucuronide", units = "umol/L", specimen = "administration site", verified = FALSE),
    jejunum_lumen_gluc  = list(analyte = "telmisartan glucuronide", units = "umol/L", specimen = "administration site", verified = FALSE),
    ileum_lumen_gluc    = list(analyte = "telmisartan glucuronide", units = "umol/L", specimen = "administration site", verified = FALSE),
    caecum_lumen_gluc   = list(analyte = "telmisartan glucuronide", units = "umol/L", specimen = "administration site", verified = FALSE),
    colon_lumen_gluc    = list(analyte = "telmisartan glucuronide", units = "umol/L", specimen = "administration site", verified = FALSE),
    duodenum_ent        = list(analyte = "telmisartan", units = "umol/L", specimen = "tissue", verified = FALSE),
    jejunum_ent         = list(analyte = "telmisartan", units = "umol/L", specimen = "tissue", verified = FALSE),
    ileum_ent           = list(analyte = "telmisartan", units = "umol/L", specimen = "tissue", verified = FALSE),
    duodenum_ent_gluc   = list(analyte = "telmisartan glucuronide", units = "umol/L", specimen = "tissue", verified = FALSE),
    jejunum_ent_gluc    = list(analyte = "telmisartan glucuronide", units = "umol/L", specimen = "tissue", verified = FALSE),
    ileum_ent_gluc      = list(analyte = "telmisartan glucuronide", units = "umol/L", specimen = "tissue", verified = FALSE),
    duodenum_muc        = list(analyte = "telmisartan", units = "umol/L", specimen = "whole blood", verified = FALSE),
    jejunum_muc         = list(analyte = "telmisartan", units = "umol/L", specimen = "whole blood", verified = FALSE),
    ileum_muc           = list(analyte = "telmisartan", units = "umol/L", specimen = "whole blood", verified = FALSE),
    serosa              = list(analyte = "telmisartan", units = "umol/L", specimen = "tissue", verified = FALSE),
    muscle              = list(analyte = "telmisartan", units = "umol/L", specimen = "tissue", verified = FALSE),
    skin                = list(analyte = "telmisartan", units = "umol/L", specimen = "tissue", verified = FALSE),
    adipose             = list(analyte = "telmisartan", units = "umol/L", specimen = "tissue", verified = FALSE),
    ehc1                = list(analyte = "telmisartan", units = "umol", specimen = "bile", verified = FALSE),
    ehc2                = list(analyte = "telmisartan", units = "umol", specimen = "bile", verified = FALSE),
    ehc3                = list(analyte = "telmisartan", units = "umol", specimen = "bile", verified = FALSE),
    ehc1_gluc           = list(analyte = "telmisartan glucuronide", units = "umol", specimen = "bile", verified = FALSE),
    ehc2_gluc           = list(analyte = "telmisartan glucuronide", units = "umol", specimen = "bile", verified = FALSE),
    ehc3_gluc           = list(analyte = "telmisartan glucuronide", units = "umol", specimen = "bile", verified = FALSE),
    a_feces             = list(analyte = "telmisartan", units = "umol", specimen = "faeces", verified = FALSE),
    a_feces_gluc        = list(analyte = "telmisartan glucuronide", units = "umol", specimen = "faeces", verified = FALSE),
    a_urine             = list(analyte = "telmisartan", units = "umol", specimen = "urine", verified = FALSE),
    auc_blood           = list(analyte = "telmisartan", units = "umol/L*h", specimen = "whole blood", verified = FALSE),
    target_central      = list(analyte = "AT1 receptor", units = "umol", specimen = "not applicable", verified = FALSE),
    target_liver1       = list(analyte = "AT1 receptor", units = "umol", specimen = "not applicable", verified = FALSE),
    target_liver2       = list(analyte = "AT1 receptor", units = "umol", specimen = "not applicable", verified = FALSE),
    target_liver3       = list(analyte = "AT1 receptor", units = "umol", specimen = "not applicable", verified = FALSE),
    target_liver4       = list(analyte = "AT1 receptor", units = "umol", specimen = "not applicable", verified = FALSE),
    target_liver5       = list(analyte = "AT1 receptor", units = "umol", specimen = "not applicable", verified = FALSE),
    target_duodenum     = list(analyte = "AT1 receptor", units = "umol", specimen = "not applicable", verified = FALSE),
    target_jejunum      = list(analyte = "AT1 receptor", units = "umol", specimen = "not applicable", verified = FALSE),
    target_ileum        = list(analyte = "AT1 receptor", units = "umol", specimen = "not applicable", verified = FALSE),
    target_muscle       = list(analyte = "AT1 receptor", units = "umol", specimen = "not applicable", verified = FALSE),
    target_skin         = list(analyte = "AT1 receptor", units = "umol", specimen = "not applicable", verified = FALSE),
    target_adipose      = list(analyte = "AT1 receptor", units = "umol", specimen = "not applicable", verified = FALSE),
    complex_central     = list(analyte = "telmisartan-AT1 complex", units = "umol", specimen = "not applicable", verified = FALSE),
    complex_liver1      = list(analyte = "telmisartan-AT1 complex", units = "umol", specimen = "not applicable", verified = FALSE),
    complex_liver2      = list(analyte = "telmisartan-AT1 complex", units = "umol", specimen = "not applicable", verified = FALSE),
    complex_liver3      = list(analyte = "telmisartan-AT1 complex", units = "umol", specimen = "not applicable", verified = FALSE),
    complex_liver4      = list(analyte = "telmisartan-AT1 complex", units = "umol", specimen = "not applicable", verified = FALSE),
    complex_liver5      = list(analyte = "telmisartan-AT1 complex", units = "umol", specimen = "not applicable", verified = FALSE),
    complex_duodenum    = list(analyte = "telmisartan-AT1 complex", units = "umol", specimen = "not applicable", verified = FALSE),
    complex_jejunum     = list(analyte = "telmisartan-AT1 complex", units = "umol", specimen = "not applicable", verified = FALSE),
    complex_ileum       = list(analyte = "telmisartan-AT1 complex", units = "umol", specimen = "not applicable", verified = FALSE),
    complex_muscle      = list(analyte = "telmisartan-AT1 complex", units = "umol", specimen = "not applicable", verified = FALSE),
    complex_skin        = list(analyte = "telmisartan-AT1 complex", units = "umol", specimen = "not applicable", verified = FALSE),
    complex_adipose     = list(analyte = "telmisartan-AT1 complex", units = "umol", specimen = "not applicable", verified = FALSE),
    occupancy_central   = list(analyte = "AT1 receptor occupancy", units = "fraction", specimen = "not applicable", verified = FALSE),
    occupancy_liver1    = list(analyte = "AT1 receptor occupancy", units = "fraction", specimen = "not applicable", verified = FALSE),
    occupancy_liver2    = list(analyte = "AT1 receptor occupancy", units = "fraction", specimen = "not applicable", verified = FALSE),
    occupancy_liver3    = list(analyte = "AT1 receptor occupancy", units = "fraction", specimen = "not applicable", verified = FALSE),
    occupancy_liver4    = list(analyte = "AT1 receptor occupancy", units = "fraction", specimen = "not applicable", verified = FALSE),
    occupancy_liver5    = list(analyte = "AT1 receptor occupancy", units = "fraction", specimen = "not applicable", verified = FALSE),
    occupancy_duodenum  = list(analyte = "AT1 receptor occupancy", units = "fraction", specimen = "not applicable", verified = FALSE),
    occupancy_jejunum   = list(analyte = "AT1 receptor occupancy", units = "fraction", specimen = "not applicable", verified = FALSE),
    occupancy_ileum     = list(analyte = "AT1 receptor occupancy", units = "fraction", specimen = "not applicable", verified = FALSE),
    occupancy_muscle    = list(analyte = "AT1 receptor occupancy", units = "fraction", specimen = "not applicable", verified = FALSE),
    occupancy_skin      = list(analyte = "AT1 receptor occupancy", units = "fraction", specimen = "not applicable", verified = FALSE),
    occupancy_adipose   = list(analyte = "AT1 receptor occupancy", units = "fraction", specimen = "not applicable", verified = FALSE)
  )

  population <- list(
    species    = "human",
    dose_range = paste(
      "single oral solution 5, 10, 20, 40, 80 and 160 mg (fasted) and",
      "single intravenous 10, 20, 40, 80 and 120 mg"
    ),
    disease_state = "healthy volunteers",
    notes = paste(
      "Tsuchitani 2024 fits individual plasma telmisartan",
      "concentration-time profiles from healthy volunteers, converted to",
      "whole-blood concentrations with the blood-to-plasma ratio Rb =",
      "0.775 (Table S6), together with pooled 72-144 h faecal excretion",
      "after 40 mg oral and intravenous dosing. The 2.5 mg oral arm was",
      "excluded because most points were below the limit of",
      "quantification. The clinical data are described as 'partly",
      "published in Stangier et al.'",
      "(J Clin Pharmacol. 2000;40:1312-1322 and",
      "J Int Med Res. 2000;28:149-167); subject counts, ages, weights and",
      "sex distribution are not reported by Tsuchitani 2024 and those two",
      "source papers are not open access and are not on disk. All",
      "physiological parameters in Tables S3-S7 are scaled to a single",
      "typical 78 kg adult."
    )
  )

  ini({
    # =================================================================
    # ESTIMATED PARAMETERS -- main-text Table 1, 'Final parameters
    # median [min-max]' column, i.e. the median over the parameter sets
    # ranked below 300 after the middle-out CGNM analysis. That is the
    # set the paper uses for the sensitivity analysis and the
    # pharmacogenomic simulation ("Parameter sets below rank 300 were
    # used for simulation"). CGNM is a fixed-effects least-squares
    # method: none of these carries an IIV term.
    #
    # UNIT CONVERSION: Table 1 reports Kd and Km,OATP1B3 in nM but the
    # ODEs of Data S1 work in umol/L (uM) throughout -- in d/dt(y1) the
    # term R_complexcentral * koff has units umol/h, so the balancing
    # term y1 * fb * R_freecentral * koff / Kd requires Kd in uM. Both
    # are therefore entered as nM / 1000. This reading is confirmed by
    # Table S10, whose twelve secondary parameters all reproduce from
    # the medians below via the auxiliary equations; e.g.
    # Vmax,OATP1B3/Km,OATP1B3 = 87,692 L/h against a published
    # 86,957 [70,211-101,277] L/h, and PSinf = 92,954 L/h against a
    # published 91,950 [74,559-106,883] L/h. Km,P-gp and Km,UGT1A3 are
    # already reported in uM and are entered unchanged.
    # =================================================================
    alpha_transit <- 0.590
    label("Scaling factor applied to every intestinal transit rate constant (unitless)")                   # Table 1 alpha, median 0.590
    beta_hep <- 0.380
    label("Fraction of hepatocyte disposition that is elimination rather than sinusoidal efflux (unitless)")  # Table 1 beta, median 0.380
    lcl_deg_feces <- log(5677)
    label("Total degradation/deconjugation clearance of Tel-GLU in the intestinal lumen (L/h)")            # Table 1 CL_degfeces, median 5677
    lcl_deg_liver <- log(7.99)
    label("Degradation/deconjugation clearance of Tel-GLU in the hepatocytes (L/h)")                       # Table 1 CL_degliver, median 7.99
    lcl_glu_ent <- log(1.44)
    label("Total Tel-GLU secretion clearance from enterocyte to small-intestinal lumen (L/h)")             # Table 1 CL_glu,ent, median 1.44
    lcl_int_all <- log(35327)
    label("Overall hepatic intrinsic clearance, unbound (L/h)")                                            # Table 1 CL_int,all, median 35,327
    lkd <- log(0.000160)
    label("Dissociation constant of unbound telmisartan for the AT1 receptor (umol/L)")                    # Table 1 Kd, median 0.160 nM = 0.000160 uM
    lkm_pgp <- log(2.97)
    label("Michaelis-Menten constant of unbound telmisartan for P-gp (umol/L)")                            # Table 1 Km,P-gp, median 2.97 uM
    lkm_ugt <- log(0.990)
    label("Michaelis-Menten constant of unbound telmisartan metabolism by UGT1A3 (umol/L)")                # Table 1 Km,UGT1A3, median 0.990 uM
    lkm_oatp1b3 <- log(0.00351)
    label("Michaelis-Menten constant of unbound telmisartan uptake by OATP1B3 (umol/L)")                   # Table 1 Km,OATP1B3, median 3.51 nM = 0.00351 uM
    lps_dif_eff <- log(230)
    label("Passive efflux clearance from hepatocyte to sinusoidal compartment (L/h)")                      # Table 1 PS_difeff, median 230
    lps_dif_ent_be <- log(73.9)
    label("Passive diffusion clearance across the enterocyte basolateral membrane (L/h)")                  # Table 1 PS_difentBE, median 73.9
    lr_total <- log(13.4)
    label("Total AT1 receptor amount in the body (umol)")                                                  # Table 1 R_total, median 13.4
    lr_dif <- log(0.0600)
    label("Ratio of passive influx clearance to OATP1B3 uptake clearance into hepatocytes (unitless)")     # Table 1 R_dif, median 0.0600
    vmax_ratio_ugt <- 0.0100
    label("Intestinal-to-hepatic expression ratio of UGT1A3 (unitless)")                                   # Table 1 VmaxtoliverUGT, median 0.0100
    f_bile <- 0.0100
    label("Fraction of hepatocyte elimination that is biliary P-gp secretion rather than UGT1A3 metabolism")  # Table 1 f_bile, median 0.0100
    lkoff <- log(1.87)
    label("Dissociation rate constant of telmisartan from the AT1 receptor (1/h)")                         # Table 1 k_off, median 1.87

    # =================================================================
    # FIXED COMPOUND-DEPENDENT CONSTANTS -- Supplementary Table S6
    # ("Drug dependent parameters used in the TMDD-PBPK model of
    # telmisartan"), except km_glu_bile which that table reports as an
    # assumed value.
    # =================================================================
    lkp_adipose <- fixed(log(0.452))
    label("Adipose-to-blood partition coefficient (unitless)")                                             # Table S6 Kp,adipose 0.452, Rodgers & Rowland
    lkp_muscle <- fixed(log(0.516))
    label("Muscle-to-blood partition coefficient (unitless)")                                              # Table S6 Kp,muscle 0.516
    lkp_skin <- fixed(log(0.335))
    label("Skin-to-blood partition coefficient (unitless)")                                                # Table S6 Kp,skin 0.335
    lkp_gut <- fixed(log(0.284))
    label("Gut-to-blood partition coefficient (unitless)")                                                 # Table S6 Kp,gut 0.284
    fu_b <- fixed(0.00645)
    label("Unbound fraction of telmisartan in blood (unitless)")                                           # Table S6 fb = fp/Rb = 0.005/0.775
    fu_h <- fixed(0.00495)
    label("Unbound fraction of telmisartan in hepatocytes (unitless)")                                     # Table S6 fh, measured in this study
    fu_h_gluc <- fixed(0.0073)
    label("Unbound fraction of Tel-GLU in hepatocytes (unitless)")                                         # Table S6 fh,glu 0.0073
    fu_gut <- fixed(0.0227)
    label("Unbound fraction of telmisartan in enterocytes (unitless)")                                     # Table S6 fgut = fp/Kp,gut = 0.005/0.284
    lkm_glu_bile <- fixed(log(100))
    label("Michaelis-Menten constant of Tel-GLU biliary secretion (umol/L)")                               # Table S6 Km,glu,bile 100 uM, assumed by the authors
    ratio_glubile_pgp <- fixed(9.13)
    label("Ratio of Tel-GLU to telmisartan biliary secretion clearance (unitless)")                        # Table S10 footnote a, from Li R 2014 sandwich-cultured hepatocytes

    # =================================================================
    # CONSTANTS SET IN THE PARAMETER BLOCK OF Data S1. Each is written
    # there as a bare assignment rather than being fitted.
    # =================================================================
    cl_r <- fixed(0)
    label("Renal clearance of telmisartan (L/h)")                                                          # Data S1 `CLr=0`
    sf_kp <- fixed(1)
    label("Scaling factor applied to every tissue-to-blood partition coefficient (unitless)")              # Data S1 `SF=1`
    gamma_ent <- fixed(1)
    label("Enterocyte basolateral-to-apical passive diffusion clearance ratio (unitless)")                 # Data S1 `Gammaent=1`
    lumen_ratio <- fixed(1)
    label("Lumen-side surface-area scaling factor for enterocyte diffusion (unitless)")                    # Data S1 `LR=1`
    apical_ratio <- fixed(1)
    label("Apical-side surface-area scaling factor for enterocyte diffusion (unitless)")                   # Data S1 `AR=1`
    fu_gut_gluc <- fixed(1)
    label("Unbound fraction of Tel-GLU in enterocytes (unitless)")                                         # Data S1 `fgutglu=1`
    vmax_ratio_pgp <- fixed(0)
    label("Intestinal-to-hepatic expression ratio of P-gp (unitless)")                                     # Data S1 `VmaxtoliverPgp = 0 #Fixed to zero considering high intestinal absorption`

    # =================================================================
    # FIXED PHYSIOLOGICAL PARAMETERS -- Supplementary Table S3 (blood
    # flows and volumes for a 78 kg adult).
    #
    # q_hepatic is the 1.4-fold-adjusted hepatic blood flow: Appendix S1
    # section 1 states that the reference value 97 L/h/78 kg gives a
    # theoretical CLh ceiling of 60.9 L/h, which the observed i.v. CLh
    # (49.5-77.1 L/h across doses) exceeds, "To solve this
    # inconsistency, Qh was adjusted to 136 L/hr". Table S3 footnote b
    # records the same 1.4-fold adjustment. The value is corroborated by
    # blood-flow mass balance in d/dt(is_liver1): the sinusoidal inflows
    # q_ha + q_serosa + q_duodenum_muc + q_jejunum_muc + q_ileum_muc =
    # 33.8 + 81.5 + 1.79 + 9.83 + 8.69 = 135.6 L/h, which matches the
    # 136 L/h outflow to 0.3% (the component flows are tabulated to
    # three significant figures).
    # =================================================================
    v_central <- fixed(5.80)
    label("Central (whole-blood) compartment volume (L)")                                                  # Table S3 Central volume 5.80
    v_muscle <- fixed(33.5)
    label("Muscle volume (L)")                                                                             # Table S3 Muscle volume 33.5
    v_skin <- fixed(8.66)
    label("Skin volume (L)")                                                                               # Table S3 Skin volume 8.66
    v_adipose <- fixed(11.1)
    label("Adipose volume (L)")                                                                            # Table S3 Adipose volume 11.1
    v_serosa <- fixed(0.548)
    label("Intestinal serosal volume (L)")                                                                 # Table S3 Serosa volume 0.548
    v_is_liver <- fixed(0.522)
    label("Total hepatic extracellular (sinusoidal) volume, split over five segments (L)")                 # Table S3 liver extracellular space 0.522
    v_int_liver <- fixed(1.36)
    label("Total hepatocellular volume, split over five segments (L)")                                     # Table S3 hepatocellular space 1.36
    v_duodenum_ent <- fixed(0.0655)
    label("Duodenal enterocyte volume (L)")                                                                # Table S3 enterocyte duodenum 0.0655
    v_jejunum_ent <- fixed(0.233)
    label("Jejunal enterocyte volume (L)")                                                                 # Table S3 enterocyte jejunum 0.233
    v_ileum_ent <- fixed(0.279)
    label("Ileal enterocyte volume (L)")                                                                   # Table S3 enterocyte ileum 0.279
    v_duodenum_muc <- fixed(0.00881)
    label("Duodenal mucosal blood volume (L)")                                                             # Table S3 mucosa duodenum 0.00881
    v_jejunum_muc <- fixed(0.0312)
    label("Jejunal mucosal blood volume (L)")                                                              # Table S3 mucosa jejunum 0.0312
    v_ileum_muc <- fixed(0.0372)
    label("Ileal mucosal blood volume (L)")                                                                # Table S3 mucosa ileum 0.0372
    v_duodenum_lumen <- fixed(0.0725)
    label("Duodenal lumen volume (L)")                                                                     # Table S3 gut lumen duodenum 0.0725
    v_jejunum_lumen <- fixed(0.258)
    label("Jejunal lumen volume (L)")                                                                      # Table S3 gut lumen jejunum 0.258
    v_ileum_lumen <- fixed(0.308)
    label("Ileal lumen volume (L)")                                                                        # Table S3 gut lumen ileum 0.308
    v_caecum_lumen <- fixed(0.0475)
    label("Caecal lumen volume (L)")                                                                       # Table S3 gut lumen caecum 0.0475
    v_colon_lumen <- fixed(0.0503)
    label("Colonic lumen volume (L)")                                                                      # Table S3 gut lumen colon 0.0503
    q_hepatic <- fixed(136)
    label("Total hepatic blood flow (L/h)")                                                                # Appendix S1 section 1: 97 L/h adjusted 1.4-fold to 136
    q_ha <- fixed(33.8)
    label("Hepatic arterial blood flow (L/h)")                                                             # Table S3 hepatic artery 33.8
    q_serosa <- fixed(81.5)
    label("Intestinal serosal blood flow (L/h)")                                                           # Table S3 serosa 81.5
    q_muscle <- fixed(50.1)
    label("Muscle blood flow (L/h)")                                                                       # Table S3 muscle 50.1
    q_skin <- fixed(20)
    label("Skin blood flow (L/h)")                                                                         # Table S3 skin 20
    q_adipose <- fixed(17.4)
    label("Adipose blood flow (L/h)")                                                                      # Table S3 adipose 17.4
    q_duodenum_muc <- fixed(1.79)
    label("Duodenal mucosal blood flow (L/h)")                                                             # Table S3 mucosa duodenum 1.79
    q_jejunum_muc <- fixed(9.83)
    label("Jejunal mucosal blood flow (L/h)")                                                              # Table S3 mucosa jejunum 9.83
    q_ileum_muc <- fixed(8.69)
    label("Ileal mucosal blood flow (L/h)")                                                                # Table S3 mucosa ileum 8.69

    # =================================================================
    # FIXED REGIONAL DISTRIBUTION FACTORS AND TRANSIT RATE CONSTANTS --
    # Supplementary Tables S4 and S5.
    # =================================================================
    f_ugt_duodenum <- fixed(0.14)
    label("Fraction of intestinal UGT1A3 in the duodenum (unitless)")                                      # Table S4 FUGT duodenum 0.14
    f_ugt_jejunum <- fixed(0.32)
    label("Fraction of intestinal UGT1A3 in the jejunum (unitless)")                                       # Table S4 FUGT jejunum 0.32
    f_ugt_ileum <- fixed(0.54)
    label("Fraction of intestinal UGT1A3 in the ileum (unitless)")                                         # Table S4 FUGT ileum 0.54
    f_dif_duodenum <- fixed(0.076)
    label("Fraction of small-intestinal surface area in the duodenum (unitless)")                          # Table S4 FDif duodenum 0.076
    f_dif_jejunum <- fixed(0.395)
    label("Fraction of small-intestinal surface area in the jejunum (unitless)")                           # Table S4 FDif jejunum 0.395
    f_dif_ileum <- fixed(0.529)
    label("Fraction of small-intestinal surface area in the ileum (unitless)")                             # Table S4 FDif ileum 0.529
    f_degfeces_duodenum <- fixed(1.00e-09)
    label("Fraction of luminal microbiota in the duodenum (unitless)")                                     # Table S4 Fdegfeces duodenum 1.00E-09
    f_degfeces_jejunum <- fixed(5.50e-08)
    label("Fraction of luminal microbiota in the jejunum (unitless)")                                      # Table S4 Fdegfeces jejunum 5.50E-08
    f_degfeces_ileum <- fixed(5.50e-05)
    label("Fraction of luminal microbiota in the ileum (unitless)")                                        # Table S4 Fdegfeces ileum 5.50E-05
    f_degfeces_caecum <- fixed(0.500)
    label("Fraction of luminal microbiota in the caecum (unitless)")                                       # Table S4 Fdegfeces caecum 0.500
    f_degfeces_colon <- fixed(0.500)
    label("Fraction of luminal microbiota in the colon (unitless)")                                        # Table S4 Fdegfeces colon 0.500
    lk_bile <- fixed(log(0.605))
    label("Transit rate constant through the three-compartment bile chain (1/h)")                          # Table S5 kbile 0.605
    lkfeces_duodenum <- fixed(log(5.12))
    label("Transit rate constant out of the duodenal lumen (1/h)")                                         # Table S5 kfeces duodenum 5.12
    lkfeces_jejunum <- fixed(log(1.47))
    label("Transit rate constant out of the jejunal lumen (1/h)")                                          # Table S5 kfeces jejunum 1.47
    lkfeces_ileum <- fixed(log(1.212))
    label("Transit rate constant out of the ileal lumen (1/h)")                                            # Table S5 kfeces ileum 1.212
    lkfeces_caecum <- fixed(log(0.477))
    label("Transit rate constant out of the caecal lumen (1/h)")                                           # Table S5 kfeces caecum 0.477
    lkfeces_colon <- fixed(log(0.159))
    label("Transit rate constant out of the colonic lumen (1/h)")                                          # Table S5 kfeces colon 0.159

    # =================================================================
    # FIXED AT1-RECEPTOR DISTRIBUTION -- Supplementary Table S7.
    #
    # The Methods state that "The fraction of receptors in the central
    # compartment was fixed as 57.2% of the total receptors in the body,
    # based on the proportion of blood volume in the organs to the total
    # systemic blood volume", and Appendix S1 states the complementary
    # tissue share is 42.8%. Both follow from the Table S7 blood-volume
    # column: the listed tissues hold 0.0376 + 0.291 + 0.587 + 0.221 +
    # 1.34 = 2.477 L of the 5.80 L total blood volume (Table S3), i.e.
    # 42.7%, leaving 57.3% in the great-vessel (central) pool.
    #
    # Table S7 then fixes "the receptor expression ratio among the
    # organs ... to the product of RNA expression level (nTPM) and organ
    # volume". Those products sum to 2814 + 227,396 + 59,155 + 974,859 +
    # 383,299 = 1,647,523, so each tissue's share of the 42.8% pool is
    # its own product divided by that sum. The five values below are
    # that arithmetic and sum to 0.428 exactly.
    # =================================================================
    at1_central <- fixed(0.572)
    label("Fraction of total AT1 receptors in the central compartment (unitless)")                         # Methods: fixed as 57.2%
    at1_liver <- fixed(0.05907382659)
    label("Fraction of total AT1 receptors in the liver, split over five segments (unitless)")             # Table S7: 0.428 * 227396 / 1647523
    at1_si <- fixed(0.0007310319795)
    label("Fraction of total AT1 receptors in the small intestine, split over three segments (unitless)")  # Table S7: 0.428 * 2814 / 1647523
    at1_muscle <- fixed(0.09957492065)
    label("Fraction of total AT1 receptors in muscle (unitless)")                                          # Table S7: 0.428 * 383299 / 1647523
    at1_skin <- fixed(0.01536751839)
    label("Fraction of total AT1 receptors in skin (unitless)")                                            # Table S7: 0.428 * 59155 / 1647523
    at1_adipose <- fixed(0.2532527024)
    label("Fraction of total AT1 receptors in adipose (unitless)")                                         # Table S7: 0.428 * 974859 / 1647523

    # =================================================================
    # RESIDUAL ERROR. Tsuchitani 2024 fits the model with CGNM, a
    # fixed-effects nonlinear-least-squares method that minimises a
    # weighted sum of squared residuals (main-text Equation 3) and
    # estimates no residual-error model at all. The term below exists
    # only so the model is a valid nlmixr2 object; it is NOT an
    # estimate from the paper and must not be read as one. Same
    # convention as Aoki_2024_bosentan_pbpk.R and
    # Tsuchitani_2026_apixaban_pbpk.R.
    # =================================================================
    propSd <- fixed(0.10)
    label("Proportional residual error placeholder (fraction)")                                            # not reported by Tsuchitani 2024; placeholder only
  })

  model({
    # -----------------------------------------------------------------
    # 1. Back-transform the log-scale parameters.
    # -----------------------------------------------------------------
    cl_deg_feces <- exp(lcl_deg_feces)
    cl_deg_liver <- exp(lcl_deg_liver)
    cl_glu_ent <- exp(lcl_glu_ent)
    cl_int_all <- exp(lcl_int_all)
    kd <- exp(lkd)
    km_pgp <- exp(lkm_pgp)
    km_ugt <- exp(lkm_ugt)
    km_oatp1b3 <- exp(lkm_oatp1b3)
    ps_dif_eff <- exp(lps_dif_eff)
    ps_dif_ent_be <- exp(lps_dif_ent_be)
    r_total <- exp(lr_total)
    r_dif <- exp(lr_dif)
    koff <- exp(lkoff)
    km_glu_bile <- exp(lkm_glu_bile)
    kp_adipose <- exp(lkp_adipose)
    kp_muscle <- exp(lkp_muscle)
    kp_skin <- exp(lkp_skin)
    kp_gut <- exp(lkp_gut)
    k_bile <- exp(lk_bile)
    kfeces_duodenum <- exp(lkfeces_duodenum)
    kfeces_jejunum <- exp(lkfeces_jejunum)
    kfeces_ileum <- exp(lkfeces_ileum)
    kfeces_caecum <- exp(lkfeces_caecum)
    kfeces_colon <- exp(lkfeces_colon)

    # -----------------------------------------------------------------
    # 2. Auxiliary functions -- Data S1, "Auxiliary functions" block,
    #    transcribed line for line. These re-express the fitted
    #    aggregate parameters (cl_int_all, beta_hep, r_dif, f_bile) as
    #    the elementary transport and metabolic capacities. The
    #    definitions of beta, CL_int,all, R_dif and f_bile are also
    #    given in the Note under main-text Table 1 and agree.
    #
    #    All twelve derived quantities are tabulated by the paper in
    #    Table S10 and reproduce here; see the vignette source-trace
    #    table for the value-by-value comparison.
    # -----------------------------------------------------------------
    ps_dif_inf <- r_dif / (1 + r_dif) * cl_int_all / beta_hep      # Data S1 PSdifinf
    gamma_hep <- ps_dif_inf / ps_dif_eff                           # Data S1 Gamma
    vmax_uptake <- km_oatp1b3 / (1 + r_dif) * cl_int_all / beta_hep  # Data S1 VmaxUptake
    vmax_pgp <- km_pgp * f_bile * cl_int_all / (1 - beta_hep) *
      r_dif / (1 + r_dif) / gamma_hep                              # Data S1 VmaxPgp
    vmax_ugt <- km_ugt * (1 - f_bile) * cl_int_all / (1 - beta_hep) *
      r_dif / (1 + r_dif) / gamma_hep                              # Data S1 VmaxUGT
    vmax_pgp_ent <- vmax_pgp * vmax_ratio_pgp                      # Data S1 VmaxPgpent
    vmax_ugt_ent <- vmax_ugt * vmax_ratio_ugt                      # Data S1 VmaxUGTent
    ps_dif_ent_eb <- ps_dif_ent_be / gamma_ent                     # Data S1 PSdifentEB
    vmax_glu_bile <- km_glu_bile * ratio_glubile_pgp * f_bile *
      cl_int_all / (1 - beta_hep) * r_dif / (1 + r_dif) / gamma_hep  # Data S1 Vmaxglubile

    # Per-segment volumes. Data S1 writes the hepatic terms as
    # `1 / (0.2 * Vi)` and `1 / (0.2 * Vh)`, i.e. the tabulated total
    # hepatic extracellular and hepatocellular volumes divided evenly
    # over the five dispersion segments.
    v_is_seg <- 0.2 * v_is_liver
    v_int_seg <- 0.2 * v_int_liver

    # Receptor amount resident in each anatomical location. Data S1
    # divides the hepatic pool over the five liver segments and the
    # small-intestinal pool over the three enterocyte segments
    # (`R_total * AT1H / 5`, `R_total * AT1SI / 3` in the RO_ equations).
    r_central <- r_total * at1_central
    r_liver_seg <- r_total * at1_liver / 5
    r_si_seg <- r_total * at1_si / 3
    r_muscle <- r_total * at1_muscle
    r_skin <- r_total * at1_skin
    r_adipose <- r_total * at1_adipose

    # -----------------------------------------------------------------
    # 3. Emergent blood concentrations leaving the perfusion-limited
    #    tissues and the serosa, and the local unbound concentrations
    #    that drive receptor binding there. Data S1 writes these
    #    inline as `y15 / Kpm / SF` etc.
    # -----------------------------------------------------------------
    cb_muscle <- muscle / kp_muscle / sf_kp
    cb_skin <- skin / kp_skin / sf_kp
    cb_adipose <- adipose / kp_adipose / sf_kp
    cb_serosa <- serosa / kp_gut / sf_kp

    # -----------------------------------------------------------------
    # 4. Receptor binding. Data S1 writes association inline as
    #    `C * fb * R_free * koff / Kd`, i.e. kon = koff / Kd applied to
    #    the unbound local concentration and the free receptor amount,
    #    and dissociation as `R_complex * koff`. Binding is fully
    #    reversible; the paper states that receptor synthesis and
    #    degradation are deliberately not modelled.
    #
    #    Data S1 references R_free<location> but does not define it,
    #    because free receptor is carried as its own state initialised
    #    to the location's receptor amount -- the same convention this
    #    laboratory's Aoki_2024_bosentan_pbpk.R uses
    #    (`d/dt(target) <- unbind - bind`, `target(0) <- rtot`). The
    #    initial conditions at the end of this block supply that, and
    #    the vignette checks the resulting conservation identity
    #    target_<loc> + complex_<loc> == r_<loc> at every time point.
    # -----------------------------------------------------------------
    kon <- koff / kd

    bind_central <- kon * fu_b * central * target_central
    bind_liver1 <- kon * fu_b * is_liver1 * target_liver1
    bind_liver2 <- kon * fu_b * is_liver2 * target_liver2
    bind_liver3 <- kon * fu_b * is_liver3 * target_liver3
    bind_liver4 <- kon * fu_b * is_liver4 * target_liver4
    bind_liver5 <- kon * fu_b * is_liver5 * target_liver5
    bind_duodenum <- kon * fu_b * duodenum_ent * target_duodenum
    bind_jejunum <- kon * fu_b * jejunum_ent * target_jejunum
    bind_ileum <- kon * fu_b * ileum_ent * target_ileum
    bind_muscle <- kon * fu_b * cb_muscle * target_muscle
    bind_skin <- kon * fu_b * cb_skin * target_skin
    bind_adipose <- kon * fu_b * cb_adipose * target_adipose

    unbind_central <- koff * complex_central
    unbind_liver1 <- koff * complex_liver1
    unbind_liver2 <- koff * complex_liver2
    unbind_liver3 <- koff * complex_liver3
    unbind_liver4 <- koff * complex_liver4
    unbind_liver5 <- koff * complex_liver5
    unbind_duodenum <- koff * complex_duodenum
    unbind_jejunum <- koff * complex_jejunum
    unbind_ileum <- koff * complex_ileum
    unbind_muscle <- koff * complex_muscle
    unbind_skin <- koff * complex_skin
    unbind_adipose <- koff * complex_adipose

    # -----------------------------------------------------------------
    # 5. Hepatic transport and metabolism per dispersion segment.
    #    Sinusoid-to-hepatocyte influx is saturable OATP1B3 uptake plus
    #    passive diffusion, both driven by the unbound sinusoidal
    #    concentration; hepatocyte disposition is passive efflux back to
    #    the sinusoid plus saturable P-gp biliary secretion and
    #    saturable UGT1A3 glucuronidation, driven by the unbound
    #    hepatocyte concentration.
    # -----------------------------------------------------------------
    uptake1 <- 0.2 * fu_b * (vmax_uptake / (km_oatp1b3 + fu_b * is_liver1) + ps_dif_inf) * is_liver1
    uptake2 <- 0.2 * fu_b * (vmax_uptake / (km_oatp1b3 + fu_b * is_liver2) + ps_dif_inf) * is_liver2
    uptake3 <- 0.2 * fu_b * (vmax_uptake / (km_oatp1b3 + fu_b * is_liver3) + ps_dif_inf) * is_liver3
    uptake4 <- 0.2 * fu_b * (vmax_uptake / (km_oatp1b3 + fu_b * is_liver4) + ps_dif_inf) * is_liver4
    uptake5 <- 0.2 * fu_b * (vmax_uptake / (km_oatp1b3 + fu_b * is_liver5) + ps_dif_inf) * is_liver5

    efflux1 <- 0.2 * fu_h * ps_dif_eff * int_liver1
    efflux2 <- 0.2 * fu_h * ps_dif_eff * int_liver2
    efflux3 <- 0.2 * fu_h * ps_dif_eff * int_liver3
    efflux4 <- 0.2 * fu_h * ps_dif_eff * int_liver4
    efflux5 <- 0.2 * fu_h * ps_dif_eff * int_liver5

    elim1 <- 0.2 * fu_h * (vmax_pgp / (km_pgp + fu_h * int_liver1) + vmax_ugt / (km_ugt + fu_h * int_liver1)) * int_liver1
    elim2 <- 0.2 * fu_h * (vmax_pgp / (km_pgp + fu_h * int_liver2) + vmax_ugt / (km_ugt + fu_h * int_liver2)) * int_liver2
    elim3 <- 0.2 * fu_h * (vmax_pgp / (km_pgp + fu_h * int_liver3) + vmax_ugt / (km_ugt + fu_h * int_liver3)) * int_liver3
    elim4 <- 0.2 * fu_h * (vmax_pgp / (km_pgp + fu_h * int_liver4) + vmax_ugt / (km_ugt + fu_h * int_liver4)) * int_liver4
    elim5 <- 0.2 * fu_h * (vmax_pgp / (km_pgp + fu_h * int_liver5) + vmax_ugt / (km_ugt + fu_h * int_liver5)) * int_liver5

    # Hepatic P-gp biliary secretion of telmisartan (feeds the bile
    # chain) and UGT1A3 glucuronidation (feeds the Tel-GLU hepatocyte
    # states).
    pgp_bile <- 0.2 * fu_h * vmax_pgp *
      (int_liver1 / (km_pgp + fu_h * int_liver1) +
         int_liver2 / (km_pgp + fu_h * int_liver2) +
         int_liver3 / (km_pgp + fu_h * int_liver3) +
         int_liver4 / (km_pgp + fu_h * int_liver4) +
         int_liver5 / (km_pgp + fu_h * int_liver5))

    ugt1 <- 0.2 * (vmax_ugt / (km_ugt + fu_h * int_liver1)) * fu_h * int_liver1
    ugt2 <- 0.2 * (vmax_ugt / (km_ugt + fu_h * int_liver2)) * fu_h * int_liver2
    ugt3 <- 0.2 * (vmax_ugt / (km_ugt + fu_h * int_liver3)) * fu_h * int_liver3
    ugt4 <- 0.2 * (vmax_ugt / (km_ugt + fu_h * int_liver4)) * fu_h * int_liver4
    ugt5 <- 0.2 * (vmax_ugt / (km_ugt + fu_h * int_liver5)) * fu_h * int_liver5

    # Tel-GLU hepatocyte disposition: deconjugation back to telmisartan
    # plus saturable biliary secretion into the Tel-GLU bile chain.
    glu_out1 <- 0.2 * (cl_deg_liver + vmax_glu_bile / (km_glu_bile + fu_h_gluc * int_liver1_gluc)) * fu_h_gluc * int_liver1_gluc
    glu_out2 <- 0.2 * (cl_deg_liver + vmax_glu_bile / (km_glu_bile + fu_h_gluc * int_liver2_gluc)) * fu_h_gluc * int_liver2_gluc
    glu_out3 <- 0.2 * (cl_deg_liver + vmax_glu_bile / (km_glu_bile + fu_h_gluc * int_liver3_gluc)) * fu_h_gluc * int_liver3_gluc
    glu_out4 <- 0.2 * (cl_deg_liver + vmax_glu_bile / (km_glu_bile + fu_h_gluc * int_liver4_gluc)) * fu_h_gluc * int_liver4_gluc
    glu_out5 <- 0.2 * (cl_deg_liver + vmax_glu_bile / (km_glu_bile + fu_h_gluc * int_liver5_gluc)) * fu_h_gluc * int_liver5_gluc

    glu_bile_in <- 0.2 * fu_h_gluc * vmax_glu_bile *
      (int_liver1_gluc / (km_glu_bile + fu_h_gluc * int_liver1_gluc) +
         int_liver2_gluc / (km_glu_bile + fu_h_gluc * int_liver2_gluc) +
         int_liver3_gluc / (km_glu_bile + fu_h_gluc * int_liver3_gluc) +
         int_liver4_gluc / (km_glu_bile + fu_h_gluc * int_liver4_gluc) +
         int_liver5_gluc / (km_glu_bile + fu_h_gluc * int_liver5_gluc))

    # -----------------------------------------------------------------
    # 6. Intestinal regional distribution factors.
    #
    #    Data S1 multiplies five processes by a regional factor:
    #    passive diffusion by Fdif_*, UGT1A3 by Fugt_*, microbial
    #    deconjugation by Fdegfeces_*, intestinal P-gp by Fpgp_*, and
    #    enterocyte-to-lumen Tel-GLU secretion by Fgluent_*. Table S4
    #    tabulates only the first three. Fpgp_* is immaterial because
    #    Data S1 fixes VmaxtoliverPgp to zero, so vmax_pgp_ent -- and
    #    with it the whole intestinal P-gp Michaelis-Menten term -- is
    #    identically zero whatever Fpgp_* may be; the Fpgp_* slot in
    #    those terms below therefore carries f_dif_* as an inert
    #    placeholder and no value is asserted for it.
    #    Fgluent_* is taken to be the
    #    surface-area vector Fdif_*, on the grounds that it gates a
    #    clearance across the enterocyte apical membrane and Table S4
    #    describes Fdif as the "ratio of diffusion clearance along the
    #    intestine"; the alternative reading (the UGT vector Fugt_*,
    #    where the conjugate is formed) is checked in the vignette and
    #    changes simulated telmisartan exposure negligibly, because
    #    CL_glu,ent is the one parameter main-text Table 1 reports as
    #    wholly unidentifiable (median 1.44 L/h over a min-max range of
    #    1.21e-05 to 5.26e+05, profile-likelihood interval [NA, NA]).
    #    Recorded in the vignette Errata.
    # -----------------------------------------------------------------
    f_gluent_duodenum <- f_dif_duodenum
    f_gluent_jejunum <- f_dif_jejunum
    f_gluent_ileum <- f_dif_ileum

    # -----------------------------------------------------------------
    # 7. ODEs -- Data S1, "ODE of TMDDPBPK model" block. The Data S1
    #    state name is given at the end of every line.
    # -----------------------------------------------------------------

    # Central (whole blood). Receives hepatic venous outflow from the
    # last sinusoidal segment, loses drug to renal clearance (zero here)
    # and to the three perfusion-limited tissues, and exchanges with the
    # central receptor pool.
    d/dt(central) <- 1 / v_central *
      (q_hepatic * (is_liver5 - central) - cl_r * central -
         q_muscle * (central - cb_muscle) -
         q_skin * (central - cb_skin) -
         q_adipose * (central - cb_adipose) +
         unbind_central - bind_central)                                     # Data S1 y1
    f(central) <- 1 / v_central

    # Liver, five sinusoidal (extracellular) segments in series. Segment
    # 1 receives the hepatic artery, the serosal drainage and the three
    # mucosal-blood drainages; segments 2-5 receive the preceding
    # segment.
    d/dt(is_liver1) <- 1 / v_is_seg *
      (q_ha * central + q_serosa * cb_serosa +
         (q_duodenum_muc * duodenum_muc + q_jejunum_muc * jejunum_muc +
            q_ileum_muc * ileum_muc) -
         q_hepatic * is_liver1 - uptake1 + efflux1 +
         unbind_liver1 - bind_liver1)                                       # Data S1 y2
    d/dt(is_liver2) <- 1 / v_is_seg *
      (q_hepatic * (is_liver1 - is_liver2) - uptake2 + efflux2 +
         unbind_liver2 - bind_liver2)                                       # Data S1 y4
    d/dt(is_liver3) <- 1 / v_is_seg *
      (q_hepatic * (is_liver2 - is_liver3) - uptake3 + efflux3 +
         unbind_liver3 - bind_liver3)                                       # Data S1 y6
    d/dt(is_liver4) <- 1 / v_is_seg *
      (q_hepatic * (is_liver3 - is_liver4) - uptake4 + efflux4 +
         unbind_liver4 - bind_liver4)                                       # Data S1 y8
    d/dt(is_liver5) <- 1 / v_is_seg *
      (q_hepatic * (is_liver4 - is_liver5) - uptake5 + efflux5 +
         unbind_liver5 - bind_liver5)                                       # Data S1 y10

    # Liver, five hepatocyte (intracellular) segments. Each gains uptake
    # from and returns efflux to its own sinusoidal segment, loses drug
    # to P-gp and UGT1A3, and regains telmisartan from hepatic
    # deconjugation of Tel-GLU.
    d/dt(int_liver1) <- 1 / v_int_seg *
      (uptake1 - efflux1 - elim1 + 0.2 * cl_deg_liver * fu_h_gluc * int_liver1_gluc)  # Data S1 y3
    d/dt(int_liver2) <- 1 / v_int_seg *
      (uptake2 - efflux2 - elim2 + 0.2 * cl_deg_liver * fu_h_gluc * int_liver2_gluc)  # Data S1 y5
    d/dt(int_liver3) <- 1 / v_int_seg *
      (uptake3 - efflux3 - elim3 + 0.2 * cl_deg_liver * fu_h_gluc * int_liver3_gluc)  # Data S1 y7
    d/dt(int_liver4) <- 1 / v_int_seg *
      (uptake4 - efflux4 - elim4 + 0.2 * cl_deg_liver * fu_h_gluc * int_liver4_gluc)  # Data S1 y9
    d/dt(int_liver5) <- 1 / v_int_seg *
      (uptake5 - efflux5 - elim5 + 0.2 * cl_deg_liver * fu_h_gluc * int_liver5_gluc)  # Data S1 y11

    # Tel-GLU in the five hepatocyte segments.
    d/dt(int_liver1_gluc) <- 1 / v_int_seg * (ugt1 - glu_out1)              # Data S1 y3G
    d/dt(int_liver2_gluc) <- 1 / v_int_seg * (ugt2 - glu_out2)              # Data S1 y5G
    d/dt(int_liver3_gluc) <- 1 / v_int_seg * (ugt3 - glu_out3)              # Data S1 y7G
    d/dt(int_liver4_gluc) <- 1 / v_int_seg * (ugt4 - glu_out4)              # Data S1 y9G
    d/dt(int_liver5_gluc) <- 1 / v_int_seg * (ugt5 - glu_out5)              # Data S1 y11G

    # Telmisartan bile chain (three transit compartments, amounts).
    d/dt(ehc1) <- pgp_bile - k_bile * ehc1                                  # Data S1 y18
    d/dt(ehc2) <- k_bile * (ehc1 - ehc2)                                    # Data S1 y19
    d/dt(ehc3) <- k_bile * (ehc2 - ehc3)                                    # Data S1 y20

    # Tel-GLU bile chain (three transit compartments, amounts).
    d/dt(ehc1_gluc) <- glu_bile_in - k_bile * ehc1_gluc                     # Data S1 y18G
    d/dt(ehc2_gluc) <- k_bile * (ehc1_gluc - ehc2_gluc)                     # Data S1 y19G
    d/dt(ehc3_gluc) <- k_bile * (ehc2_gluc - ehc3_gluc)                     # Data S1 y20G

    # Small-intestinal lumen, three segments. Each receives bile (the
    # duodenum only), transit from the preceding segment, P-gp efflux
    # and passive diffusion out of the enterocyte, and telmisartan
    # regenerated by microbial deconjugation of luminal Tel-GLU; each
    # loses drug to transit and to enterocyte uptake. The oral dose
    # enters here.
    d/dt(duodenum_lumen) <- 1 / v_duodenum_lumen *
      (k_bile * ehc3 -
         (alpha_transit * kfeces_duodenum * v_duodenum_lumen +
            lumen_ratio * ps_dif_ent_be * f_dif_duodenum) * duodenum_lumen +
         fu_gut * (vmax_pgp_ent * f_dif_duodenum / (km_pgp + fu_gut * duodenum_ent) +
                     apical_ratio * ps_dif_ent_eb * f_dif_duodenum) * duodenum_ent +
         duodenum_lumen_gluc * cl_deg_feces * f_degfeces_duodenum)          # Data S1 y12
    f(duodenum_lumen) <- 1 / v_duodenum_lumen

    d/dt(jejunum_lumen) <- 1 / v_jejunum_lumen *
      (alpha_transit * kfeces_duodenum * v_duodenum_lumen * duodenum_lumen -
         (alpha_transit * kfeces_jejunum * v_jejunum_lumen +
            lumen_ratio * ps_dif_ent_be * f_dif_jejunum) * jejunum_lumen +
         fu_gut * (vmax_pgp_ent * f_dif_jejunum / (km_pgp + fu_gut * jejunum_ent) +
                     apical_ratio * ps_dif_ent_eb * f_dif_jejunum) * jejunum_ent +
         jejunum_lumen_gluc * cl_deg_feces * f_degfeces_jejunum)            # Data S1 y13

    d/dt(ileum_lumen) <- 1 / v_ileum_lumen *
      (alpha_transit * kfeces_jejunum * v_jejunum_lumen * jejunum_lumen -
         (alpha_transit * kfeces_ileum * v_ileum_lumen +
            lumen_ratio * ps_dif_ent_be * f_dif_ileum) * ileum_lumen +
         fu_gut * (vmax_pgp_ent * f_dif_ileum / (km_pgp + fu_gut * ileum_ent) +
                     apical_ratio * ps_dif_ent_eb * f_dif_ileum) * ileum_ent +
         ileum_lumen_gluc * cl_deg_feces * f_degfeces_ileum)                # Data S1 y14

    # Tel-GLU in the small-intestinal lumen. Receives Tel-GLU bile and
    # enterocyte secretion, loses it to transit and microbial
    # deconjugation.
    d/dt(duodenum_lumen_gluc) <- 1 / v_duodenum_lumen *
      (k_bile * ehc3_gluc -
         alpha_transit * kfeces_duodenum * v_duodenum_lumen * duodenum_lumen_gluc -
         duodenum_lumen_gluc * cl_deg_feces * f_degfeces_duodenum +
         fu_gut_gluc * duodenum_ent_gluc * cl_glu_ent * f_gluent_duodenum)  # Data S1 y12G

    d/dt(jejunum_lumen_gluc) <- 1 / v_jejunum_lumen *
      (alpha_transit * kfeces_duodenum * v_duodenum_lumen * duodenum_lumen_gluc -
         alpha_transit * kfeces_jejunum * v_jejunum_lumen * jejunum_lumen_gluc -
         jejunum_lumen_gluc * cl_deg_feces * f_degfeces_jejunum +
         fu_gut_gluc * jejunum_ent_gluc * cl_glu_ent * f_gluent_jejunum)    # Data S1 y13G

    d/dt(ileum_lumen_gluc) <- 1 / v_ileum_lumen *
      (alpha_transit * kfeces_jejunum * v_jejunum_lumen * jejunum_lumen_gluc -
         ileum_lumen_gluc * cl_deg_feces * f_degfeces_ileum -
         alpha_transit * kfeces_ileum * v_ileum_lumen * ileum_lumen_gluc +
         fu_gut_gluc * ileum_ent_gluc * cl_glu_ent * f_gluent_ileum)        # Data S1 y14G

    # Enterocytes, three segments. Take up drug from the lumen, lose it
    # back to the lumen by P-gp and passive diffusion, exchange with the
    # mucosal blood, metabolise it by intestinal UGT1A3, and bind the
    # enterocyte AT1 receptor pool.
    d/dt(duodenum_ent) <- 1 / v_duodenum_ent *
      (lumen_ratio * ps_dif_ent_be * f_dif_duodenum * duodenum_lumen -
         fu_gut * (vmax_pgp_ent * f_dif_duodenum / (km_pgp + fu_gut * duodenum_ent) +
                     apical_ratio * ps_dif_ent_eb * f_dif_duodenum) * duodenum_ent -
         fu_gut * ps_dif_ent_eb * f_dif_duodenum * duodenum_ent +
         fu_b * ps_dif_ent_be * f_dif_duodenum * duodenum_muc -
         fu_gut * vmax_ugt_ent * f_ugt_duodenum / (km_ugt + fu_gut * duodenum_ent) * duodenum_ent +
         unbind_duodenum - bind_duodenum)                                   # Data S1 y21

    d/dt(jejunum_ent) <- 1 / v_jejunum_ent *
      (lumen_ratio * ps_dif_ent_be * f_dif_jejunum * jejunum_lumen -
         fu_gut * (vmax_pgp_ent * f_dif_jejunum / (km_pgp + fu_gut * jejunum_ent) +
                     apical_ratio * ps_dif_ent_eb * f_dif_jejunum) * jejunum_ent -
         fu_gut * ps_dif_ent_eb * f_dif_jejunum * jejunum_ent +
         fu_b * ps_dif_ent_be * f_dif_jejunum * jejunum_muc -
         fu_gut * vmax_ugt_ent * f_ugt_jejunum / (km_ugt + fu_gut * jejunum_ent) * jejunum_ent +
         unbind_jejunum - bind_jejunum)                                     # Data S1 y23

    d/dt(ileum_ent) <- 1 / v_ileum_ent *
      (lumen_ratio * ps_dif_ent_be * f_dif_ileum * ileum_lumen -
         fu_gut * (vmax_pgp_ent * f_dif_ileum / (km_pgp + fu_gut * ileum_ent) +
                     apical_ratio * ps_dif_ent_eb * f_dif_ileum) * ileum_ent -
         fu_gut * ps_dif_ent_eb * f_dif_ileum * ileum_ent +
         fu_b * ps_dif_ent_be * f_dif_ileum * ileum_muc -
         fu_gut * vmax_ugt_ent * f_ugt_ileum / (km_ugt + fu_gut * ileum_ent) * ileum_ent +
         unbind_ileum - bind_ileum)                                         # Data S1 y25

    # Tel-GLU in the enterocytes: formed by intestinal UGT1A3, secreted
    # to the lumen.
    d/dt(duodenum_ent_gluc) <- 1 / v_duodenum_ent *
      (fu_gut * duodenum_ent * vmax_ugt_ent * f_ugt_duodenum / (km_ugt + fu_gut * duodenum_ent) -
         fu_gut_gluc * duodenum_ent_gluc * cl_glu_ent * f_gluent_duodenum)  # Data S1 y21G
    d/dt(jejunum_ent_gluc) <- 1 / v_jejunum_ent *
      (fu_gut * jejunum_ent * vmax_ugt_ent * f_ugt_jejunum / (km_ugt + fu_gut * jejunum_ent) -
         fu_gut_gluc * jejunum_ent_gluc * cl_glu_ent * f_gluent_jejunum)    # Data S1 y23G
    d/dt(ileum_ent_gluc) <- 1 / v_ileum_ent *
      (fu_gut * ileum_ent * vmax_ugt_ent * f_ugt_ileum / (km_ugt + fu_gut * ileum_ent) -
         fu_gut_gluc * ileum_ent_gluc * cl_glu_ent * f_gluent_ileum)        # Data S1 y25G

    # Mucosal blood, three segments. Perfused from the central
    # compartment, drains to the liver sinusoid, exchanges with its
    # enterocyte by passive diffusion.
    d/dt(duodenum_muc) <- 1 / v_duodenum_muc *
      (q_duodenum_muc * (central - duodenum_muc) +
         fu_gut * ps_dif_ent_eb * f_dif_duodenum * duodenum_ent -
         fu_b * ps_dif_ent_be * f_dif_duodenum * duodenum_muc)              # Data S1 y22
    d/dt(jejunum_muc) <- 1 / v_jejunum_muc *
      (q_jejunum_muc * (central - jejunum_muc) +
         fu_gut * ps_dif_ent_eb * f_dif_jejunum * jejunum_ent -
         fu_b * ps_dif_ent_be * f_dif_jejunum * jejunum_muc)                # Data S1 y24
    d/dt(ileum_muc) <- 1 / v_ileum_muc *
      (q_ileum_muc * (central - ileum_muc) +
         fu_gut * ps_dif_ent_eb * f_dif_ileum * ileum_ent -
         fu_b * ps_dif_ent_be * f_dif_ileum * ileum_muc)                    # Data S1 y26

    # Serosa (the non-absorptive limb of the segregated-flow intestine).
    d/dt(serosa) <- 1 / v_serosa * (q_serosa * (central - cb_serosa))       # Data S1 y27

    # Perfusion-limited tissues, each carrying its own AT1 receptor pool.
    d/dt(muscle) <- 1 / v_muscle *
      (q_muscle * (central - cb_muscle) + unbind_muscle - bind_muscle)      # Data S1 y15
    d/dt(skin) <- 1 / v_skin *
      (q_skin * (central - cb_skin) + unbind_skin - bind_skin)              # Data S1 y16
    d/dt(adipose) <- 1 / v_adipose *
      (q_adipose * (central - cb_adipose) + unbind_adipose - bind_adipose)  # Data S1 y17

    # Caecum and colon: transit only, with microbial interconversion
    # between telmisartan and Tel-GLU, then irreversible faecal output.
    d/dt(caecum_lumen) <- 1 / v_caecum_lumen *
      (alpha_transit * kfeces_ileum * v_ileum_lumen * ileum_lumen +
         caecum_lumen_gluc * cl_deg_feces * f_degfeces_caecum -
         alpha_transit * kfeces_caecum * v_caecum_lumen * caecum_lumen)     # Data S1 y29
    d/dt(colon_lumen) <- 1 / v_colon_lumen *
      (alpha_transit * kfeces_caecum * v_caecum_lumen * caecum_lumen +
         colon_lumen_gluc * cl_deg_feces * f_degfeces_colon -
         alpha_transit * colon_lumen * v_colon_lumen * kfeces_colon)        # Data S1 y34
    d/dt(a_feces) <- alpha_transit * colon_lumen * v_colon_lumen * kfeces_colon  # Data S1 y35

    d/dt(caecum_lumen_gluc) <- 1 / v_caecum_lumen *
      (alpha_transit * kfeces_ileum * v_ileum_lumen * ileum_lumen_gluc -
         alpha_transit * kfeces_caecum * v_caecum_lumen * caecum_lumen_gluc -
         caecum_lumen_gluc * cl_deg_feces * f_degfeces_caecum)              # Data S1 y29G
    d/dt(colon_lumen_gluc) <- 1 / v_colon_lumen *
      (alpha_transit * kfeces_caecum * v_caecum_lumen * caecum_lumen_gluc -
         colon_lumen_gluc * cl_deg_feces * f_degfeces_colon -
         alpha_transit * colon_lumen_gluc * v_colon_lumen * kfeces_colon)   # Data S1 y34G
    d/dt(a_feces_gluc) <- alpha_transit * colon_lumen_gluc * v_colon_lumen * kfeces_colon  # Data S1 y35G

    # Cumulative blood AUC and urinary output (the latter is
    # identically zero because cl_r is fixed to zero).
    d/dt(auc_blood) <- central                                              # Data S1 y31
    d/dt(a_urine) <- central * cl_r                                         # Data S1 y33

    # Free AT1 receptor in each location. Data S1 uses R_free* without
    # defining it; it is the complement of the bound complex within a
    # conserved pool, initialised below.
    d/dt(target_central) <- unbind_central - bind_central
    d/dt(target_liver1) <- unbind_liver1 - bind_liver1
    d/dt(target_liver2) <- unbind_liver2 - bind_liver2
    d/dt(target_liver3) <- unbind_liver3 - bind_liver3
    d/dt(target_liver4) <- unbind_liver4 - bind_liver4
    d/dt(target_liver5) <- unbind_liver5 - bind_liver5
    d/dt(target_duodenum) <- unbind_duodenum - bind_duodenum
    d/dt(target_jejunum) <- unbind_jejunum - bind_jejunum
    d/dt(target_ileum) <- unbind_ileum - bind_ileum
    d/dt(target_muscle) <- unbind_muscle - bind_muscle
    d/dt(target_skin) <- unbind_skin - bind_skin
    d/dt(target_adipose) <- unbind_adipose - bind_adipose

    # Telmisartan-AT1 complex in each location.
    d/dt(complex_central) <- bind_central - unbind_central                  # Data S1 R_complexcentral
    d/dt(complex_liver1) <- bind_liver1 - unbind_liver1                     # Data S1 R_complexhep1
    d/dt(complex_liver2) <- bind_liver2 - unbind_liver2                     # Data S1 R_complexhep2
    d/dt(complex_liver3) <- bind_liver3 - unbind_liver3                     # Data S1 R_complexhep3
    d/dt(complex_liver4) <- bind_liver4 - unbind_liver4                     # Data S1 R_complexhep4
    d/dt(complex_liver5) <- bind_liver5 - unbind_liver5                     # Data S1 R_complexhep5
    d/dt(complex_duodenum) <- bind_duodenum - unbind_duodenum               # Data S1 R_complexduodenum
    d/dt(complex_jejunum) <- bind_jejunum - unbind_jejunum                  # Data S1 R_complexjejunum
    d/dt(complex_ileum) <- bind_ileum - unbind_ileum                        # Data S1 R_complexileum
    d/dt(complex_muscle) <- bind_muscle - unbind_muscle                     # Data S1 R_complexmuscle
    d/dt(complex_skin) <- bind_skin - unbind_skin                           # Data S1 R_complexskin
    d/dt(complex_adipose) <- bind_adipose - unbind_adipose                  # Data S1 R_complexadipose

    # Fractional AT1 receptor occupancy in each location, carried as an
    # explicit state exactly as Data S1 carries it.
    d/dt(occupancy_central) <- (bind_central - unbind_central) / r_central        # Data S1 RO_central
    d/dt(occupancy_liver1) <- (bind_liver1 - unbind_liver1) / r_liver_seg         # Data S1 RO_hep1
    d/dt(occupancy_liver2) <- (bind_liver2 - unbind_liver2) / r_liver_seg         # Data S1 RO_hep2
    d/dt(occupancy_liver3) <- (bind_liver3 - unbind_liver3) / r_liver_seg         # Data S1 RO_hep3
    d/dt(occupancy_liver4) <- (bind_liver4 - unbind_liver4) / r_liver_seg         # Data S1 RO_hep4
    d/dt(occupancy_liver5) <- (bind_liver5 - unbind_liver5) / r_liver_seg         # Data S1 RO_hep5
    d/dt(occupancy_duodenum) <- (bind_duodenum - unbind_duodenum) / r_si_seg      # Data S1 RO_mucosa1
    d/dt(occupancy_jejunum) <- (bind_jejunum - unbind_jejunum) / r_si_seg         # Data S1 RO_mucosa2
    d/dt(occupancy_ileum) <- (bind_ileum - unbind_ileum) / r_si_seg               # Data S1 RO_mucosa3
    d/dt(occupancy_muscle) <- (bind_muscle - unbind_muscle) / r_muscle            # Data S1 RO_muscle
    d/dt(occupancy_skin) <- (bind_skin - unbind_skin) / r_skin                    # Data S1 RO_skin
    d/dt(occupancy_adipose) <- (bind_adipose - unbind_adipose) / r_adipose        # Data S1 RO_adipose

    # -----------------------------------------------------------------
    # 8. Initial conditions. Every drug state starts empty; every
    #    receptor pool starts fully unbound at its location's share of
    #    r_total.
    # -----------------------------------------------------------------
    target_central(0) <- r_central
    target_liver1(0) <- r_liver_seg
    target_liver2(0) <- r_liver_seg
    target_liver3(0) <- r_liver_seg
    target_liver4(0) <- r_liver_seg
    target_liver5(0) <- r_liver_seg
    target_duodenum(0) <- r_si_seg
    target_jejunum(0) <- r_si_seg
    target_ileum(0) <- r_si_seg
    target_muscle(0) <- r_muscle
    target_skin(0) <- r_skin
    target_adipose(0) <- r_adipose

    # -----------------------------------------------------------------
    # 9. Observation. The paper converts observed plasma concentrations
    #    to whole-blood concentrations with Rb before fitting, so the
    #    central compartment is whole blood and Cc is a whole-blood
    #    concentration in umol/L.
    # -----------------------------------------------------------------
    Cc <- central
    Cc ~ prop(propSd)
  })
}
