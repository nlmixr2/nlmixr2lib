Tsuchitani_2026_apixaban_pbpk <- function() {
  description <- paste(
    "PBPK (whole-body, 47 ODEs, top-down CGNM fit). Apixaban disposition",
    "after single intravenous (0.5-5 mg) and oral (0.5-50 mg) doses in",
    "healthy adults, built by Tsuchitani et al. (2026) to test whether",
    "biliary secretion and enterohepatic circulation (EHC) - long believed",
    "negligible for apixaban - in fact drive its faecal elimination. Four",
    "mechanistic blocks are coupled: (1) a segregated-flow intestinal model",
    "in which duodenum, jejunum and ileum each carry a luminal, an",
    "enterocyte and a mucosal-blood space, with the serosal blood returning",
    "separately to the liver; (2) a membrane-limited five-compartment",
    "tandem-dispersion liver in which each hepatic extracellular",
    "(sinusoidal) sub-compartment exchanges with its own hepatocyte by",
    "passive diffusion only (apixaban is not an OATP1B substrate, so the",
    "saturable uptake term is present but its Vmax is fixed to zero);",
    "(3) a mechanistic kidney resolving glomerulus, proximal tubule,",
    "distal tubule and collecting duct, each as a vascular / epithelial-cell",
    "/ tubular-lumen triple; and (4) a three-compartment EHC transit chain",
    "carrying canalicular bile through the biliary tree and gallbladder",
    "back into the duodenal lumen at rate k_bile. Hepatic elimination is",
    "split by f_bile between P-gp-mediated biliary secretion and CYP3A4",
    "metabolism; intestinal P-gp and CYP3A4 are scaled from the hepatic",
    "Vmax values. Eleven parameters were estimated by the Cluster",
    "Gauss-Newton method against blood concentrations plus urinary, faecal",
    "and total-metabolite mass-balance data; five (CL_int_all, f_bile,",
    "PS_dif_ent_BE, SF and k_bile) were identifiable by profile likelihood.",
    "Values below are the medians of the 200 bootstrap parameter sets,",
    "which is the set the paper's own Discussion quotes. CGNM is a",
    "fixed-effects method, so the model carries no between-subject",
    "variability and no residual-error estimate; propSd is a placeholder.",
    "The activated-charcoal interaction sub-model of the same paper is NOT",
    "included - see the vignette Errata."
  )
  reference <- paste(
    "Tsuchitani T, Kou W, Tomi M, Sugiyama Y. Characterizing apixaban",
    "pharmacokinetics through physiologically-based pharmacokinetic",
    "modeling: critical role of biliary secretion and enterohepatic",
    "circulation in humans. CPT Pharmacometrics Syst Pharmacol.",
    "2026;15:e70163. doi:10.1002/psp4.70163. (Open access, CC BY-NC;",
    "accepted 12 November 2025, published in volume 15.)",
    "The ODE system and every fixed constant are transcribed from",
    "Supporting Information Data S1 (psp470163-sup-0001-supinfo.txt), the",
    "authors' model code; the physiological constants are independently",
    "tabulated in Data S2 (psp470163-sup-0002-supinfo.pdf) Tables S3A,",
    "S3B, S3C and S4, and every one of them agrees with the code. The",
    "eleven estimated parameters are the 'Median' column of main-text",
    "Table 1."
  )
  vignette <- "Tsuchitani_2026_apixaban_pbpk"

  # Time in hours; the published code works in molar units throughout.
  # Tissue states hold concentrations in umol/L (micromolar) and the
  # luminal / excretion / bile states hold amounts in umol, so doses must
  # be supplied in umol. The paper does not print apixaban's molecular
  # weight, but its own Table S1A pins it: the 5 mg iv row gives
  # AUCinf = 3.43 uM*h and Dose/AUCinf = 3.17 L/h, so the dose is
  # 3.43 * 3.17 = 10.87 umol and MW = 5000 / 10.87 = 460 g/mol (the
  # 2.5 and 1.25 mg rows give 458 and 458). See the vignette.
  units <- list(time = "h", dosing = "umol", concentration = "umol/L")

  # Canonical states used here: `depot` (oral dose in the stomach),
  # `central` (blood), `is_liver1..5` / `int_liver1..5` (the hepatic
  # extracellular and hepatocyte sub-compartments of
  # `pbpkSubCompartmentRegex`, exactly as in Aoki_2024_bosentan_pbpk.R
  # from the same laboratory), `muscle`, `skin`, `adipose`, and
  # `a_urine`.
  #
  # Everything below is genuinely paper-mechanistic. The segregated-flow
  # intestine needs THREE spaces per gut segment (lumen / enterocyte /
  # mucosal blood) plus a shared serosal space, which no canonical
  # covers - the register's bare `duodenum` / `jejunum` / `ileum` are
  # single lumped segments. The mechanistic kidney needs a
  # vascular / cell / tubule triple for each of four nephron segments.
  # The EHC transit chain is a three-compartment delay, not the
  # register's gated `gallbladder`. Following the precedent set by
  # Aoki_2024_bosentan_pbpk.R, these are declared paper-specific rather
  # than silently registered as ~25 new canonicals.
  paper_specific_compartments <- c(
    "depot_iv",
    "is_liver1", "is_liver2", "is_liver3", "is_liver4", "is_liver5",
    "int_liver1", "int_liver2", "int_liver3", "int_liver4", "int_liver5",
    "duodenum_lumen", "jejunum_lumen", "ileum_lumen",
    "duodenum_ent", "jejunum_ent", "ileum_ent",
    "duodenum_muc", "jejunum_muc", "ileum_muc",
    "serosa",
    "ehc1", "ehc2", "ehc3",
    "blood_gl", "tube_gl",
    "blood_pt", "cell_pt", "tube_pt",
    "blood_dt", "cell_dt", "tube_dt",
    "blood_cd", "cell_cd", "tube_cd",
    "a_feces", "a_metab_gut", "a_metab_liver",
    "a_bile_canalicular", "a_bile_duodenum",
    "auc_blood"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. verified = FALSE means NOT checked against the
  # source paper.
  compartmentData <- list(
    depot              = list(analyte = "apixaban", units = "umol", specimen = "administration site", verified = FALSE),
    depot_iv           = list(analyte = "apixaban", units = "umol", specimen = "administration site", verified = FALSE),
    central            = list(analyte = "apixaban", units = "umol/L", specimen = "whole blood", verified = FALSE),
    is_liver1          = list(analyte = "apixaban", units = "umol/L", specimen = "tissue", verified = FALSE),
    is_liver2          = list(analyte = "apixaban", units = "umol/L", specimen = "tissue", verified = FALSE),
    is_liver3          = list(analyte = "apixaban", units = "umol/L", specimen = "tissue", verified = FALSE),
    is_liver4          = list(analyte = "apixaban", units = "umol/L", specimen = "tissue", verified = FALSE),
    is_liver5          = list(analyte = "apixaban", units = "umol/L", specimen = "tissue", verified = FALSE),
    int_liver1         = list(analyte = "apixaban", units = "umol/L", specimen = "tissue", verified = FALSE),
    int_liver2         = list(analyte = "apixaban", units = "umol/L", specimen = "tissue", verified = FALSE),
    int_liver3         = list(analyte = "apixaban", units = "umol/L", specimen = "tissue", verified = FALSE),
    int_liver4         = list(analyte = "apixaban", units = "umol/L", specimen = "tissue", verified = FALSE),
    int_liver5         = list(analyte = "apixaban", units = "umol/L", specimen = "tissue", verified = FALSE),
    duodenum_lumen     = list(analyte = "apixaban", units = "umol/L", specimen = "administration site", verified = FALSE),
    jejunum_lumen      = list(analyte = "apixaban", units = "umol/L", specimen = "administration site", verified = FALSE),
    ileum_lumen        = list(analyte = "apixaban", units = "umol/L", specimen = "administration site", verified = FALSE),
    duodenum_ent       = list(analyte = "apixaban", units = "umol/L", specimen = "tissue", verified = FALSE),
    jejunum_ent        = list(analyte = "apixaban", units = "umol/L", specimen = "tissue", verified = FALSE),
    ileum_ent          = list(analyte = "apixaban", units = "umol/L", specimen = "tissue", verified = FALSE),
    duodenum_muc       = list(analyte = "apixaban", units = "umol/L", specimen = "whole blood", verified = FALSE),
    jejunum_muc        = list(analyte = "apixaban", units = "umol/L", specimen = "whole blood", verified = FALSE),
    ileum_muc          = list(analyte = "apixaban", units = "umol/L", specimen = "whole blood", verified = FALSE),
    serosa             = list(analyte = "apixaban", units = "umol/L", specimen = "tissue", verified = FALSE),
    muscle             = list(analyte = "apixaban", units = "umol/L", specimen = "tissue", verified = FALSE),
    skin               = list(analyte = "apixaban", units = "umol/L", specimen = "tissue", verified = FALSE),
    adipose            = list(analyte = "apixaban", units = "umol/L", specimen = "tissue", verified = FALSE),
    ehc1               = list(analyte = "apixaban", units = "umol", specimen = "bile", verified = FALSE),
    ehc2               = list(analyte = "apixaban", units = "umol", specimen = "bile", verified = FALSE),
    ehc3               = list(analyte = "apixaban", units = "umol", specimen = "bile", verified = FALSE),
    blood_gl           = list(analyte = "apixaban", units = "umol/L", specimen = "whole blood", verified = FALSE),
    tube_gl            = list(analyte = "apixaban", units = "umol/L", specimen = "urine", verified = FALSE),
    blood_pt           = list(analyte = "apixaban", units = "umol/L", specimen = "whole blood", verified = FALSE),
    cell_pt            = list(analyte = "apixaban", units = "umol/L", specimen = "tissue", verified = FALSE),
    tube_pt            = list(analyte = "apixaban", units = "umol/L", specimen = "urine", verified = FALSE),
    blood_dt           = list(analyte = "apixaban", units = "umol/L", specimen = "whole blood", verified = FALSE),
    cell_dt            = list(analyte = "apixaban", units = "umol/L", specimen = "tissue", verified = FALSE),
    tube_dt            = list(analyte = "apixaban", units = "umol/L", specimen = "urine", verified = FALSE),
    blood_cd           = list(analyte = "apixaban", units = "umol/L", specimen = "whole blood", verified = FALSE),
    cell_cd            = list(analyte = "apixaban", units = "umol/L", specimen = "tissue", verified = FALSE),
    tube_cd            = list(analyte = "apixaban", units = "umol/L", specimen = "urine", verified = FALSE),
    a_urine            = list(analyte = "apixaban", units = "umol", specimen = "urine", verified = FALSE),
    a_feces            = list(analyte = "apixaban", units = "umol", specimen = "faeces", verified = FALSE),
    a_metab_gut        = list(analyte = "apixaban metabolites", units = "umol", specimen = "not applicable", verified = FALSE),
    a_metab_liver      = list(analyte = "apixaban metabolites", units = "umol", specimen = "not applicable", verified = FALSE),
    a_bile_canalicular = list(analyte = "apixaban", units = "umol", specimen = "bile", verified = FALSE),
    a_bile_duodenum    = list(analyte = "apixaban", units = "umol", specimen = "bile", verified = FALSE),
    auc_blood          = list(analyte = "apixaban", units = "umol/L*h", specimen = "whole blood", verified = FALSE)
  )

  # No covariates. Every volume, flow and partition coefficient in Data S1
  # is a fixed constant for a single 78 kg reference adult; the published
  # code carries no body-weight or demographic scaling of any kind.
  covariateData <- list()

  population <- list(
    species    = "human",
    weight     = "78 kg reference adult (the single body weight hard-coded as BW in Data S1; every volume and flow in Tables S3A and S3B is reported per 78 kg)",
    dose_range = "0.5, 1.25, 2.5, 3.75 and 5 mg single intravenous doses; 0.5, 1, 2.5, 5, 10, 20, 25 and 50 mg single oral doses",
    notes      = paste(
      "The model was fitted top-down to digitised published data rather",
      "than to individual-level records: blood concentrations after single",
      "iv and oral doses were extracted from Frost et al. and Raghavan",
      "et al. with WebPlotDigitizer 4.7 and converted from plasma to",
      "total blood with Rb = 1.09 (Methods 2.2). Urinary excretion after",
      "each iv dose and after 20 mg po, plus faecal excretion and total",
      "generated metabolites after 20 mg po, were used as additional",
      "fitting targets, with a weight of 3 on the faecal and",
      "total-metabolite series. Data beyond 60 h post-dose were excluded",
      "because sampling schedules differed across dose levels. Subject",
      "counts, ages, sexes and weights of the underlying healthy-volunteer",
      "studies are not restated by Tsuchitani 2026 and belong to the",
      "primary publications."
    )
  )

  ini({
    # =================================================================
    # ESTIMATED PARAMETERS -- main-text Table 1, 'Median' column of the
    # 200 bootstrap parameter sets. The 'Best fit' column is NOT usable
    # as a complete set: it reports NA for six of the eleven parameters
    # (only a profile-likelihood interval is given for those), whereas
    # the bootstrap median column is complete. The paper's own
    # Discussion quotes the median column ("CL_int,all (22.5 L/h,
    # Table 1)" and "f_bile (0.44, Table 1)"), so this is the reference
    # parameter set of the publication. CGNM is a fixed-effects
    # least-squares method: none of these carries an IIV term.
    # =================================================================
    beta_hep <- 0.983
    label("Fraction of hepatocyte disposition that is elimination rather than sinusoidal efflux (unitless)")      # Table 1 beta, median 0.983
    lcl_int_all <- log(22.5)
    label("Overall hepatic intrinsic clearance, unbound (L/h)")                                                   # Table 1 CL_int,all, median 22.5
    lcl_pgp_renal <- log(1.29)
    label("P-gp-mediated transport clearance in the renal proximal tubule (L/h)")                                 # Table 1 CL_Pgp,renal, median 1.29
    f_bile <- 0.441
    label("Fraction of hepatic intrinsic clearance that is biliary secretion rather than CYP3A4 metabolism")      # Table 1 f_bile, median 0.441
    lk_bile <- log(0.241)
    label("Transit rate constant through the three-compartment bile chain (1/h)")                                 # Table 1 k_bile, median 0.241
    llr_renal <- log(6.41)
    label("Apical-to-basolateral membrane surface-area ratio in proximal tubule epithelium (unitless)")           # Table 1 LR_renal, median 6.41
    lp_d <- log(3.4e-07)
    label("Passive diffusion clearance per unit membrane surface area in renal cells (L/h/m^2)")                  # Table 1 P_d, median 3.4e-07
    lps_dif_ent_be <- log(5.88)
    label("Passive diffusion clearance across the enterocyte basolateral membrane (L/h)")                         # Table 1 PS_difentBE, median 5.88
    lsf_kp <- log(0.235)
    label("Scaling factor applied to every tissue-to-blood partition coefficient (unitless)")                     # Table 1 SF, median 0.235
    vmax_ratio_cyp <- 0.00165
    label("Intestinal-to-hepatic expression ratio of CYP3A4 (unitless)")                                          # Table 1 VmaxtoliverCYP, median 0.00165
    vmax_ratio_pgp <- 0.183
    label("Intestinal-to-hepatic expression ratio of P-gp (unitless)")                                            # Table 1 VmaxtoliverPgp, median 0.183

    # =================================================================
    # FIXED COMPOUND-DEPENDENT CONSTANTS -- Table S4 and the parameter
    # block of Data S1. Every value below appears in both sources and
    # they agree.
    # =================================================================
    km_cyp <- fixed(100)
    label("Michaelis constant for CYP3A4 metabolism of unbound apixaban (umol/L)")                                # Table S4 KmCYP = 100 (Wang 2010); Data S1 KmCYP=100
    km_pgp <- fixed(100)
    label("Michaelis constant for P-gp transport of unbound apixaban (umol/L)")                                   # Table S4 KmPgp = 100 (Zhang 2013); Data S1 KmPgp=100
    km_uptake <- fixed(100)
    label("Michaelis constant for hepatic basolateral active uptake (umol/L)")                                    # Data S1 KmUptake=100 (inactive: vmax_uptake is 0)
    vmax_uptake <- fixed(0)
    label("Maximum velocity of hepatic basolateral active uptake (umol/h)")                                       # Data S1 VmaxUptake = 0; apixaban is not an OATP1B1/3 substrate (Methods 2.1)
    rb <- fixed(1.09)
    label("Blood-to-plasma concentration ratio (unitless)")                                                       # Table S4 Rb = 1.09; Data S1 Rb=1.09
    fu_p <- fixed(0.131)
    label("Unbound fraction in plasma (unitless)")                                                                # Table S4 fp = 0.131; Data S1 fp=0.131
    fu_h <- fixed(0.239)
    label("Unbound fraction in hepatocytes (unitless)")                                                           # Table S4 fh = 0.239; Data S1 fh=0.239
    fu_gut <- fixed(0.065)
    label("Unbound fraction in enterocytes (unitless)")                                                           # Table S4 fgut = 0.065; Data S1 fgut=0.065
    fu_r <- fixed(0.115)
    label("Unbound fraction in renal epithelial cells (unitless)")                                                # Table S4 fr = 0.115; Data S1 fr=0.115
    lkp_adipose <- fixed(log(1.372))
    label("Log adipose-to-blood partition coefficient before SF scaling (unitless)")                              # Data S1 Kpa=1.372 (Table S4 prints 1.37; Rodgers & Rowland 2006)
    lkp_muscle <- fixed(log(0.363))
    label("Log muscle-to-blood partition coefficient before SF scaling (unitless)")                               # Data S1 Kpm=0.363; Table S4 Kpm = 0.363
    lkp_skin <- fixed(log(1.325))
    label("Log skin-to-blood partition coefficient before SF scaling (unitless)")                                 # Data S1 Kps=1.325 (Table S4 prints 1.32)
    lkp_gut <- fixed(log(0.980))
    label("Log serosa-to-blood partition coefficient before SF scaling (unitless)")                               # Data S1 Kpgut=0.980; Table S4 Kpgut = 0.98
    gamma_ent <- fixed(1)
    label("Enterocyte diffusion influx-to-efflux clearance ratio (unitless)")                                     # Table S4 Gammaent = 1, assumed (apixaban is neutral)
    lr_ent <- fixed(1)
    label("Apical-to-basolateral ratio of enterocyte diffusion influx clearance (unitless)")                      # Table S4 LR = 1, assumed; Data S1 LR=1
    ar_ent <- fixed(1)
    label("Apical-to-basolateral ratio of enterocyte diffusion efflux clearance (unitless)")                      # Table S4 AR = 1, assumed; Data S1 AR=LR
    lka <- fixed(log(4))
    label("Gastric emptying rate constant, stomach to duodenal lumen (1/h)")                                      # Table S3C ka = 4 /h (Jamei 2009); Data S1 ka=4
    k_dose <- fixed(100)
    label("First-order rate constant approximating the intravenous bolus (1/h)")                                  # Data S1 kdose=100

    # =================================================================
    # FIXED PHYSIOLOGICAL CONSTANTS -- Tables S3A / S3B / S3C. Data S1
    # writes each as a coefficient times BW = 78 kg; the products below
    # reproduce the Table S3 values to the printed precision (verified
    # value by value in the vignette source-trace table).
    # =================================================================
    # --- Volumes (L) -------------------------------------------------
    v_central <- fixed(5.7954)
    label("Central (blood) volume (L)")                                                                           # Data S1 Vcentral=0.0743*78; Table S3A 5.8
    v_adipose <- fixed(11.076)
    label("Adipose volume (L)")                                                                                   # Data S1 Va=142*78/1000; Table S3A 11.1
    v_muscle <- fixed(33.462)
    label("Muscle volume (L)")                                                                                    # Data S1 Vm=429*78/1000; Table S3A 33.5
    v_skin <- fixed(8.658)
    label("Skin volume (L)")                                                                                      # Data S1 Vs=111*78/1000; Table S3A 8.66
    v_serosa <- fixed(0.54756)
    label("Intestinal serosa volume (L)")                                                                         # Data S1 Vsero=0.54756; Table S3A 0.548 (Yoshikado 2021)
    v_liver_cell <- fixed(1.3572)
    label("Total hepatocellular volume (L)")                                                                      # Data S1 Vh=17.4*78/1000; Table S3A 1.36
    v_liver_ex <- fixed(0.52182)
    label("Total hepatic extracellular (sinusoidal) volume (L)")                                                  # Data S1 Vi=6.69*78/1000; Table S3A 0.522
    v_duodenum_lumen <- fixed(0.07254)
    label("Duodenal luminal volume (L)")                                                                          # Data S1 Vduodenum_lumen=0.07254; Table S3A 0.0725 (GastroPlus)
    v_jejunum_lumen <- fixed(0.25818)
    label("Jejunal luminal volume (L)")                                                                           # Data S1 Vjejunum_lumen=0.25818; Table S3A 0.258
    v_ileum_lumen <- fixed(0.3081)
    label("Ileal luminal volume (L)")                                                                             # Data S1 Vileum_lumen=0.3081; Table S3A 0.308
    v_duodenum_ent <- fixed(0.06552)
    label("Duodenal enterocyte volume (L)")                                                                       # Data S1 Vduodenum_ent=0.06552; Table S3A 0.0655 (Simcyp 19.1)
    v_jejunum_ent <- fixed(0.23322)
    label("Jejunal enterocyte volume (L)")                                                                        # Data S1 Vjejunum_ent=0.23322; Table S3A 0.233
    v_ileum_ent <- fixed(0.27846)
    label("Ileal enterocyte volume (L)")                                                                          # Data S1 Vileum_ent=0.27846; Table S3A 0.278
    v_duodenum_muc <- fixed(0.008814)
    label("Duodenal mucosal blood volume (L)")                                                                    # Data S1 Vduodenum_muc=0.008814; Table S3A 0.00881 (Kawai 1998)
    v_jejunum_muc <- fixed(0.0312)
    label("Jejunal mucosal blood volume (L)")                                                                     # Data S1 Vjejunum_muc=0.0312; Table S3A 0.0312
    v_ileum_muc <- fixed(0.037206)
    label("Ileal mucosal blood volume (L)")                                                                       # Data S1 Vileum_muc=0.037206; Table S3A 0.0372
    v_blood_gl <- fixed(0.003627)
    label("Glomerular vascular volume (L)")                                                                       # Data S1 VbloodGL=0.0000465*78; Table S3A 0.00363 (Nishiyama 2019)
    v_tube_gl <- fixed(0.0055614)
    label("Glomerular tubular (filtrate) volume (L)")                                                             # Data S1 VtubeGL=0.0000713*78; Table S3A 0.00556
    v_blood_pt <- fixed(0.036972)
    label("Proximal-tubule vascular volume (L)")                                                                  # Data S1 VbloodPT=0.000474*78; Table S3A 0.037
    v_cell_pt <- fixed(0.03081)
    label("Proximal-tubule epithelial-cell volume (L)")                                                           # Data S1 VcellPT=0.000395*78; Table S3A 0.0308
    v_tube_pt <- fixed(0.056706)
    label("Proximal-tubule luminal volume (L)")                                                                   # Data S1 VtubePT=0.000727*78; Table S3A 0.0567
    v_blood_dt <- fixed(0.007956)
    label("Distal-tubule vascular volume (L)")                                                                    # Data S1 VbloodDT=0.000102*78; Table S3A 0.00796
    v_cell_dt <- fixed(0.0066456)
    label("Distal-tubule epithelial-cell volume (L)")                                                             # Data S1 VcellDT=0.0000852*78; Table S3A 0.00665
    v_tube_dt <- fixed(0.012246)
    label("Distal-tubule luminal volume (L)")                                                                     # Data S1 VtubeDT=0.000157*78; Table S3A 0.0122
    v_blood_cd <- fixed(0.045006)
    label("Collecting-duct vascular volume (L)")                                                                  # Data S1 VbloodCD=0.000577*78; Table S3A 0.045
    v_cell_cd <- fixed(0.03744)
    label("Collecting-duct epithelial-cell volume (L)")                                                           # Data S1 VcellCD=0.000480*78; Table S3A 0.0374
    v_tube_cd <- fixed(0.068952)
    label("Collecting-duct luminal volume (L)")                                                                   # Data S1 VtubeCD=0.000884*78; Table S3A 0.069

    # --- Blood and urine flows (L/h) ---------------------------------
    q_adipose <- fixed(17.4096)
    label("Adipose blood flow (L/h)")                                                                             # Data S1 Qa=3.72*78*60/1000; Table S3B 17.4
    q_muscle <- fixed(50.076)
    label("Muscle blood flow (L/h)")                                                                              # Data S1 Qm=10.7*78*60/1000; Table S3B 50.1
    q_skin <- fixed(20.0304)
    label("Skin blood flow (L/h)")                                                                                # Data S1 Qs=4.28*78*60/1000; Table S3B 20
    q_liver <- fixed(96.72)
    label("Total hepatic blood flow (L/h)")                                                                       # Data S1 Qh=(0.0164+0.09+0.0796+0.744+0.310)*78; Table S3B 96.7
    q_hepatic_artery <- fixed(24.18)
    label("Hepatic arterial blood flow (L/h)")                                                                    # Data S1 Qha=0.310*78; Table S3B 24.2
    q_serosa <- fixed(58.032)
    label("Intestinal serosal blood flow (L/h)")                                                                  # Data S1 Qsero=0.744*78; Table S3B 58
    q_duodenum_muc <- fixed(1.2792)
    label("Duodenal mucosal blood flow (L/h)")                                                                    # Data S1 Qduodenum_muc=0.0164*78; Table S3B 1.28 (Pang 2020)
    q_jejunum_muc <- fixed(7.02)
    label("Jejunal mucosal blood flow (L/h)")                                                                     # Data S1 Qjejunum_muc=0.09*78; Table S3B 7.02
    q_ileum_muc <- fixed(6.2088)
    label("Ileal mucosal blood flow (L/h)")                                                                       # Data S1 Qileum_muc=0.0796*78; Table S3B 6.21
    q_renal <- fixed(82.68)
    label("Total renal blood flow (L/h)")                                                                         # Data S1 Qr=1.06*78; Table S3B 82.7
    q_gfr <- fixed(8.346)
    label("Glomerular filtration rate (L/h)")                                                                     # Data S1 QGFR=0.107*78; Table S3B 8.35
    q_urine_pt_dt <- fixed(2.886)
    label("Tubular fluid flow, proximal to distal tubule (L/h)")                                                  # Data S1 QuPTtoDT=0.0370*78; Table S3B 2.89 (Scotcher 2016)
    q_urine_dt_cd <- fixed(0.77532)
    label("Tubular fluid flow, distal tubule to collecting duct (L/h)")                                           # Data S1 QuDTtoCD=0.00994*78; Table S3B 0.775
    q_urine_cd <- fixed(0.064974)
    label("Urine flow out of the collecting duct (L/h)")                                                          # Data S1 QuCD=0.000833*78; Table S3B 0.065

    # --- Renal membrane surface areas (m^2) --------------------------
    sa_v_pt <- fixed(0.81)
    label("Basolateral membrane surface area, proximal tubule (m^2)")                                             # Table S3C SAvpt = 0.81 (Scotcher 2016); Data S1 SAvpt=0.81
    sa_u_pt <- fixed(6.1)
    label("Apical membrane surface area, proximal tubule (m^2)")                                                  # Table S3C SAupt = 6.1; Data S1 SAupt=6.1
    sa_v_dt <- fixed(0.21)
    label("Basolateral membrane surface area, distal tubule (m^2)")                                               # Table S3C SAvdt = 0.21; Data S1 SAvdt=0.21
    sa_u_dt <- fixed(0.21)
    label("Apical membrane surface area, distal tubule (m^2)")                                                    # Table S3C SAudt = 0.21; Data S1 SAudt=0.21
    sa_v_cd <- fixed(0.045)
    label("Basolateral membrane surface area, collecting duct (m^2)")                                             # Table S3C SAvcd = 0.045; Data S1 SAvcd=0.045
    sa_u_cd <- fixed(0.045)
    label("Apical membrane surface area, collecting duct (m^2)")                                                  # Table S3C SAucd = 0.045; Data S1 SAucd=0.045

    # --- Regional enzyme / transporter / permeability fractions ------
    f_cyp_duodenum <- fixed(0.14)
    label("Fraction of intestinal CYP3A4 in the duodenum (unitless)")                                             # Table S3C Fcyp_duodenum = 0.14 (Simcyp 19.1)
    f_cyp_jejunum <- fixed(0.32)
    label("Fraction of intestinal CYP3A4 in the jejunum (unitless)")                                              # Table S3C Fcyp_jejunum = 0.32
    f_cyp_ileum <- fixed(0.54)
    label("Fraction of intestinal CYP3A4 in the ileum (unitless)")                                                # Table S3C Fcyp_ileum = 0.54
    f_pgp_duodenum <- fixed(0.0197)
    label("Fraction of intestinal P-gp in the duodenum (unitless)")                                               # Table S3C Fpgp_duodenum = 0.0197
    f_pgp_jejunum <- fixed(0.371)
    label("Fraction of intestinal P-gp in the jejunum (unitless)")                                                # Table S3C Fpgp_jejunum = 0.371
    f_pgp_ileum <- fixed(0.61)
    label("Fraction of intestinal P-gp in the ileum (unitless)")                                                  # Table S3C Fpgp_ileum = 0.61
    f_dif_duodenum <- fixed(0.076)
    label("Fraction of intestinal passive diffusion clearance in the duodenum (unitless)")                        # Table S3C Fdif_duodenum = 0.076
    f_dif_jejunum <- fixed(0.395)
    label("Fraction of intestinal passive diffusion clearance in the jejunum (unitless)")                         # Table S3C Fdif_jejunum = 0.395
    f_dif_ileum <- fixed(0.529)
    label("Fraction of intestinal passive diffusion clearance in the ileum (unitless)")                           # Table S3C Fdif_ileum = 0.529

    # --- Intestinal transit ------------------------------------------
    alpha_transit <- fixed(0.5)
    label("Scaling factor applied to every intestinal transit rate constant (unitless)")                          # Table S3C Alpha = 0.5; Data S1 Alpha=0.5
    k_feces_duodenum <- fixed(5.12)
    label("Transit rate constant, duodenum to jejunum, before alpha scaling (1/h)")                               # Data S1 kfeces_duodenum=2.56*2; Table S3C 5.12 (GastroPlus)
    k_feces_jejunum <- fixed(1.212)
    label("Transit rate constant, jejunum to ileum, before alpha scaling (1/h)")                                  # Data S1 kfeces_jejunum=0.606*2; Table S3C 1.21
    k_feces_ileum <- fixed(1.47)
    label("Transit rate constant, ileum to faeces, before alpha scaling (1/h)")                                   # Data S1 kfeces_ileum=0.735*2; Table S3C 1.47

    # =================================================================
    # RESIDUAL ERROR
    # CGNM minimises a weighted sum of squared residuals (Equation 1)
    # and estimates neither between-subject variability nor a residual
    # variance; the paper reports no error model at all. nlmixr2 model
    # definitions require a residual-error term, so propSd below is a
    # fixed placeholder for syntactic completeness only and must NOT be
    # read as an estimate. Same convention as
    # Aoki_2024_bosentan_pbpk.R, Mi_2023_cefquinome_pbpk.R and
    # An_2012_mitoxantrone_human_pbpk.R.
    # =================================================================
    propSd <- fixed(0.10)
    label("Proportional residual error placeholder (fraction)")                                                   # not reported by Tsuchitani 2026; placeholder only
  })

  model({
    # -----------------------------------------------------------------
    # 1. Back-transform the log-scale parameters.
    # -----------------------------------------------------------------
    cl_int_all <- exp(lcl_int_all)
    cl_pgp_renal <- exp(lcl_pgp_renal)
    k_bile <- exp(lk_bile)
    lr_renal <- exp(llr_renal)
    p_d <- exp(lp_d)
    ps_dif_ent_be <- exp(lps_dif_ent_be)
    sf_kp <- exp(lsf_kp)
    kp_adipose <- exp(lkp_adipose)
    kp_muscle <- exp(lkp_muscle)
    kp_skin <- exp(lkp_skin)
    kp_gut <- exp(lkp_gut)
    ka <- exp(lka)

    # -----------------------------------------------------------------
    # 2. Derived compound constants -- the header block of Data S1.
    #
    #    The hepatic parameterisation is defined by the note under
    #    Table 1: beta = (CL_Pgp + CL_CYP) / (PS_difeff + CL_Pgp +
    #    CL_CYP), CL_int,all = PS_difinf * beta, and f_bile = CL_Pgp /
    #    (CL_Pgp + CL_CYP). Inverting those three relations gives the
    #    three assignments below, which is exactly what Data S1 codes.
    # -----------------------------------------------------------------
    fu_b <- fu_p / rb                                     # Data S1 fb=fp/Rb; Table S4 fb = 0.12
    ps_dif_inf <- cl_int_all / beta_hep                   # Data S1 PSdifinf = CLintall / Beta
    ps_dif_eff <- cl_int_all / beta_hep                   # Data S1 PSdifeff = CLintall / Beta
    vmax_pgp <- km_pgp * f_bile / (1 - beta_hep) * cl_int_all         # Data S1 VmaxPgp
    vmax_cyp <- km_cyp * (1 - f_bile) / (1 - beta_hep) * cl_int_all   # Data S1 VmaxCYP
    vmax_pgp_ent <- vmax_pgp * vmax_ratio_pgp             # Data S1 VmaxPgpent = VmaxPgp * VmaxtoliverPgp
    vmax_cyp_ent <- vmax_cyp * vmax_ratio_cyp             # Data S1 VmaxCYPent = VmaxCYP * VmaxtoliverCYP
    ps_dif_ent_eb <- ps_dif_ent_be / gamma_ent            # Data S1 PSdifentEB = PSdifentBE / Gammaent

    #    The liver is split into five equal tandem sub-compartments, so
    #    each holds a fraction fdisp = 1/5 of the hepatic volume and
    #    receives fdisp of every liver-wide clearance. Data S1 writes
    #    this fraction as the literal 0.2 throughout the hepatic ODEs.
    fdisp <- 0.2

    #    Renal permeability-surface-area products. With the fitted
    #    P_d of 3.4e-07 L/h/m^2 every one of these is numerically
    #    negligible, which is the paper's finding that apixaban
    #    undergoes neither tubular secretion nor reabsorption and that
    #    CL_renal is simply fu_b * GFR (Discussion; Table S1A).
    ps_dif_pro_ce_ve <- p_d * sa_v_pt                     # Data S1 PSdifproCEtoVE
    ps_dif_pro_ve_ce <- p_d * sa_v_pt                     # Data S1 PSdifproVEtoCE
    ps_dif_pro_ce_pt <- p_d * sa_u_pt * lr_renal          # Data S1 PSdifproCEtoPT
    ps_dif_pro_pt_ce <- p_d * sa_u_pt * lr_renal          # Data S1 PSdifproPTtoCE
    ps_dif_dis_ce_ve <- p_d * sa_v_dt                     # Data S1 PSdifdisCEtoVE
    ps_dif_dis_ve_ce <- p_d * sa_v_dt                     # Data S1 PSdifdisVEtoCE
    ps_dif_dis_ce_dt <- p_d * sa_u_dt                     # Data S1 PSdifdisCEtoDT
    ps_dif_dis_dt_ce <- p_d * sa_u_dt                     # Data S1 PSdifdisDTtoCE
    ps_dif_col_ce_ve <- p_d * sa_v_cd                     # Data S1 PSdifcolCEtoVE
    ps_dif_col_ve_ce <- p_d * sa_v_cd                     # Data S1 PSdifcolVEtoCE
    ps_dif_col_ce_cd <- p_d * sa_u_cd                     # Data S1 PSdifcolCEtoCD
    ps_dif_col_cd_ce <- p_d * sa_u_cd                     # Data S1 PSdifcolCDtoCE

    #    Renal plasma flows between nephron segments.
    q_renal_gl_pt <- q_renal - q_gfr                                          # Data S1 QrGLtoPT
    q_renal_pt_dt <- q_renal - 2 * q_gfr / 3 - q_urine_cd / 3                 # Data S1 QrPTtoDT
    q_renal_dt_cd <- q_renal - q_gfr / 3 - 2 * q_urine_cd / 3                 # Data S1 QrDTtoCD

    # -----------------------------------------------------------------
    # 3. Repeated flux expressions. Naming them once keeps the ODEs
    #    below a line-for-line match to Data S1.
    # -----------------------------------------------------------------
    #    Canalicular P-gp secretion out of each hepatocyte
    #    sub-compartment, summed over the five (umol/h).
    bile_flux <- fdisp * fu_h * vmax_pgp *
      (int_liver1 / (km_pgp + fu_h * int_liver1) +
         int_liver2 / (km_pgp + fu_h * int_liver2) +
         int_liver3 / (km_pgp + fu_h * int_liver3) +
         int_liver4 / (km_pgp + fu_h * int_liver4) +
         int_liver5 / (km_pgp + fu_h * int_liver5))

    #    Sinusoid-to-hepatocyte influx (saturable uptake + passive) and
    #    hepatocyte-to-sinusoid passive efflux, per sub-compartment.
    upt1 <- fu_b * (vmax_uptake / (km_uptake + fu_b * is_liver1) + ps_dif_inf) * is_liver1
    upt2 <- fu_b * (vmax_uptake / (km_uptake + fu_b * is_liver2) + ps_dif_inf) * is_liver2
    upt3 <- fu_b * (vmax_uptake / (km_uptake + fu_b * is_liver3) + ps_dif_inf) * is_liver3
    upt4 <- fu_b * (vmax_uptake / (km_uptake + fu_b * is_liver4) + ps_dif_inf) * is_liver4
    upt5 <- fu_b * (vmax_uptake / (km_uptake + fu_b * is_liver5) + ps_dif_inf) * is_liver5
    eff1 <- fu_h * ps_dif_eff * int_liver1
    eff2 <- fu_h * ps_dif_eff * int_liver2
    eff3 <- fu_h * ps_dif_eff * int_liver3
    eff4 <- fu_h * ps_dif_eff * int_liver4
    eff5 <- fu_h * ps_dif_eff * int_liver5

    #    Total intracellular elimination (P-gp to bile + CYP3A4) out of
    #    each hepatocyte sub-compartment.
    elim1 <- fu_h * (vmax_pgp / (km_pgp + fu_h * int_liver1) + vmax_cyp / (km_cyp + fu_h * int_liver1)) * int_liver1
    elim2 <- fu_h * (vmax_pgp / (km_pgp + fu_h * int_liver2) + vmax_cyp / (km_cyp + fu_h * int_liver2)) * int_liver2
    elim3 <- fu_h * (vmax_pgp / (km_pgp + fu_h * int_liver3) + vmax_cyp / (km_cyp + fu_h * int_liver3)) * int_liver3
    elim4 <- fu_h * (vmax_pgp / (km_pgp + fu_h * int_liver4) + vmax_cyp / (km_cyp + fu_h * int_liver4)) * int_liver4
    elim5 <- fu_h * (vmax_pgp / (km_pgp + fu_h * int_liver5) + vmax_cyp / (km_cyp + fu_h * int_liver5)) * int_liver5

    #    Hepatic CYP3A4 metabolite formation, summed over the five
    #    hepatocyte sub-compartments (umol/h).
    met_liver_flux <- fdisp * fu_h *
      (vmax_cyp / (km_cyp + fu_h * int_liver1) * int_liver1 +
         vmax_cyp / (km_cyp + fu_h * int_liver2) * int_liver2 +
         vmax_cyp / (km_cyp + fu_h * int_liver3) * int_liver3 +
         vmax_cyp / (km_cyp + fu_h * int_liver4) * int_liver4 +
         vmax_cyp / (km_cyp + fu_h * int_liver5) * int_liver5)

    #    Enterocyte CYP3A4 metabolite formation, per gut segment.
    met_duo <- fu_gut * vmax_cyp_ent * f_cyp_duodenum / (km_cyp + fu_gut * duodenum_ent) * duodenum_ent
    met_jej <- fu_gut * vmax_cyp_ent * f_cyp_jejunum / (km_cyp + fu_gut * jejunum_ent) * jejunum_ent
    met_ile <- fu_gut * vmax_cyp_ent * f_cyp_ileum / (km_cyp + fu_gut * ileum_ent) * ileum_ent

    #    Enterocyte-to-lumen efflux (apical P-gp + apical passive), per
    #    gut segment (umol/h).
    sec_duo <- fu_gut * (vmax_pgp_ent * f_pgp_duodenum / (km_pgp + fu_gut * duodenum_ent) +
                           ar_ent * ps_dif_ent_eb * f_dif_duodenum) * duodenum_ent
    sec_jej <- fu_gut * (vmax_pgp_ent * f_pgp_jejunum / (km_pgp + fu_gut * jejunum_ent) +
                           ar_ent * ps_dif_ent_eb * f_dif_jejunum) * jejunum_ent
    sec_ile <- fu_gut * (vmax_pgp_ent * f_pgp_ileum / (km_pgp + fu_gut * ileum_ent) +
                           ar_ent * ps_dif_ent_eb * f_dif_ileum) * ileum_ent

    #    Delivery of bile into the duodenal lumen (umol/h).
    bile_to_duodenum <- k_bile * ehc3

    # -----------------------------------------------------------------
    # 4. ODE system -- transcribed line for line from Data S1. Tissue
    #    states hold CONCENTRATIONS (umol/L); the dose, bile, faecal,
    #    urinary and metabolite states hold AMOUNTS (umol), exactly as
    #    the published code does.
    # -----------------------------------------------------------------

    #    Dosing compartments. The oral dose empties from the stomach
    #    into the duodenal lumen at ka; the intravenous dose enters
    #    blood at k_dose = 100 /h, a numerical stand-in for a bolus
    #    (half-life 25 s).
    d/dt(depot) <- -depot * ka                            # Data S1 d/dt(XpoDose)
    d/dt(depot_iv) <- -depot_iv * k_dose                  # Data S1 d/dt(Xdose)

    #    Blood. Receives hepatic venous return (is_liver5) and renal
    #    venous return (blood_cd); loses total hepatic and renal flow
    #    and exchanges with the three perfusion-limited tissues.
    d/dt(central) <- (depot_iv * k_dose +
      q_liver * (is_liver5 - central) -
      q_renal * central + q_renal * blood_cd -
      q_muscle * (central - muscle / kp_muscle / sf_kp) -
      q_skin * (central - skin / kp_skin / sf_kp) -
      q_adipose * (central - adipose / kp_adipose / sf_kp)) / v_central   # Data S1 d/dt(y1)

    #    Hepatic extracellular (sinusoidal) tandem chain. The first
    #    sub-compartment receives the hepatic artery plus all splanchnic
    #    venous return (serosa and the three mucosal beds); the rest
    #    receive their predecessor's outflow.
    d/dt(is_liver1) <- (q_hepatic_artery * central +
      q_serosa * serosa / kp_gut / sf_kp +
      q_duodenum_muc * duodenum_muc + q_jejunum_muc * jejunum_muc + q_ileum_muc * ileum_muc -
      q_liver * is_liver1 - fdisp * upt1 + fdisp * eff1) / (fdisp * v_liver_ex)   # Data S1 d/dt(y2)
    d/dt(is_liver2) <- (q_liver * (is_liver1 - is_liver2) -
      fdisp * upt2 + fdisp * eff2) / (fdisp * v_liver_ex)                         # Data S1 d/dt(y4)
    d/dt(is_liver3) <- (q_liver * (is_liver2 - is_liver3) -
      fdisp * upt3 + fdisp * eff3) / (fdisp * v_liver_ex)                         # Data S1 d/dt(y6)
    d/dt(is_liver4) <- (q_liver * (is_liver3 - is_liver4) -
      fdisp * upt4 + fdisp * eff4) / (fdisp * v_liver_ex)                         # Data S1 d/dt(y8)
    d/dt(is_liver5) <- (q_liver * (is_liver4 - is_liver5) -
      fdisp * upt5 + fdisp * eff5) / (fdisp * v_liver_ex)                         # Data S1 d/dt(y10)

    #    Hepatocyte chain: passive uptake in, passive efflux out, and
    #    intracellular elimination by canalicular P-gp and CYP3A4.
    d/dt(int_liver1) <- (fdisp * upt1 - fdisp * eff1 - fdisp * elim1) / (fdisp * v_liver_cell)   # Data S1 d/dt(y3)
    d/dt(int_liver2) <- (fdisp * upt2 - fdisp * eff2 - fdisp * elim2) / (fdisp * v_liver_cell)   # Data S1 d/dt(y5)
    d/dt(int_liver3) <- (fdisp * upt3 - fdisp * eff3 - fdisp * elim3) / (fdisp * v_liver_cell)   # Data S1 d/dt(y7)
    d/dt(int_liver4) <- (fdisp * upt4 - fdisp * eff4 - fdisp * elim4) / (fdisp * v_liver_cell)   # Data S1 d/dt(y9)
    d/dt(int_liver5) <- (fdisp * upt5 - fdisp * eff5 - fdisp * elim5) / (fdisp * v_liver_cell)   # Data S1 d/dt(y11)

    #    Enterohepatic circulation: a three-compartment transit chain
    #    from the canaliculus (ehc1) through the biliary tree and
    #    gallbladder (ehc2, ehc3) into the duodenal lumen. Amounts.
    d/dt(ehc1) <- bile_flux - k_bile * ehc1                # Data S1 d/dt(y18)
    d/dt(ehc2) <- k_bile * (ehc1 - ehc2)                   # Data S1 d/dt(y19)
    d/dt(ehc3) <- k_bile * (ehc2 - ehc3)                   # Data S1 d/dt(y20)

    #    Intestinal lumen (segregated-flow model). Each segment receives
    #    transit from upstream, loses drug to downstream transit and to
    #    basolateral-direction absorption into its enterocyte, and gains
    #    apically secreted drug back from that enterocyte. The duodenum
    #    additionally receives gastric emptying and bile.
    d/dt(duodenum_lumen) <- (depot * ka + bile_to_duodenum -
      (alpha_transit * k_feces_duodenum * v_duodenum_lumen +
         lr_ent * ps_dif_ent_be * f_dif_duodenum) * duodenum_lumen +
      sec_duo) / v_duodenum_lumen                          # Data S1 d/dt(y12)
    d/dt(jejunum_lumen) <- (alpha_transit * k_feces_duodenum * v_duodenum_lumen * duodenum_lumen -
      (alpha_transit * k_feces_jejunum * v_jejunum_lumen +
         lr_ent * ps_dif_ent_be * f_dif_jejunum) * jejunum_lumen +
      sec_jej) / v_jejunum_lumen                           # Data S1 d/dt(y13)
    d/dt(ileum_lumen) <- (alpha_transit * k_feces_jejunum * v_jejunum_lumen * jejunum_lumen -
      (alpha_transit * k_feces_ileum * v_ileum_lumen +
         lr_ent * ps_dif_ent_be * f_dif_ileum) * ileum_lumen +
      sec_ile) / v_ileum_lumen                             # Data S1 d/dt(y14)

    #    Enterocytes: absorption in from the lumen, apical secretion and
    #    CYP3A4 metabolism out, and passive exchange with mucosal blood.
    d/dt(duodenum_ent) <- (lr_ent * ps_dif_ent_be * f_dif_duodenum * duodenum_lumen - sec_duo -
      fu_gut * ps_dif_ent_eb * f_dif_duodenum * duodenum_ent +
      fu_b * ps_dif_ent_be * f_dif_duodenum * duodenum_muc - met_duo) / v_duodenum_ent   # Data S1 d/dt(y21)
    d/dt(jejunum_ent) <- (lr_ent * ps_dif_ent_be * f_dif_jejunum * jejunum_lumen - sec_jej -
      fu_gut * ps_dif_ent_eb * f_dif_jejunum * jejunum_ent +
      fu_b * ps_dif_ent_be * f_dif_jejunum * jejunum_muc - met_jej) / v_jejunum_ent      # Data S1 d/dt(y23)
    d/dt(ileum_ent) <- (lr_ent * ps_dif_ent_be * f_dif_ileum * ileum_lumen - sec_ile -
      fu_gut * ps_dif_ent_eb * f_dif_ileum * ileum_ent +
      fu_b * ps_dif_ent_be * f_dif_ileum * ileum_muc - met_ile) / v_ileum_ent            # Data S1 d/dt(y25)

    #    Mucosal blood beds: perfused from blood, drained to the liver.
    d/dt(duodenum_muc) <- (q_duodenum_muc * (central - duodenum_muc) +
      fu_gut * ps_dif_ent_eb * f_dif_duodenum * duodenum_ent -
      fu_b * ps_dif_ent_be * f_dif_duodenum * duodenum_muc) / v_duodenum_muc   # Data S1 d/dt(y22)
    d/dt(jejunum_muc) <- (q_jejunum_muc * (central - jejunum_muc) +
      fu_gut * ps_dif_ent_eb * f_dif_jejunum * jejunum_ent -
      fu_b * ps_dif_ent_be * f_dif_jejunum * jejunum_muc) / v_jejunum_muc      # Data S1 d/dt(y24)
    d/dt(ileum_muc) <- (q_ileum_muc * (central - ileum_muc) +
      fu_gut * ps_dif_ent_eb * f_dif_ileum * ileum_ent -
      fu_b * ps_dif_ent_be * f_dif_ileum * ileum_muc) / v_ileum_muc            # Data S1 d/dt(y26)

    #    Intestinal serosa: a perfusion-limited splanchnic tissue that
    #    drains to the liver alongside the mucosal beds.
    d/dt(serosa) <- q_serosa * (central - serosa / kp_gut / sf_kp) / v_serosa  # Data S1 d/dt(y27)

    #    Perfusion-limited, well-stirred peripheral tissues.
    d/dt(muscle) <- q_muscle * (central - muscle / kp_muscle / sf_kp) / v_muscle      # Data S1 d/dt(y15)
    d/dt(skin) <- q_skin * (central - skin / kp_skin / sf_kp) / v_skin                # Data S1 d/dt(y16)
    d/dt(adipose) <- q_adipose * (central - adipose / kp_adipose / sf_kp) / v_adipose # Data S1 d/dt(y17)

    #    Mechanistic kidney (Nishiyama 2019 / Scotcher 2016 topology).
    #    Blood is filtered at the glomerulus; the filtrate passes down
    #    the tubular lumen while the epithelial cells of each segment
    #    exchange passively with both the vascular and the luminal side.
    #    Apical P-gp in the proximal tubule is the only active term.
    d/dt(blood_gl) <- (q_renal * central - q_renal_gl_pt * blood_gl -
      fu_b * q_gfr * blood_gl) / v_blood_gl                                    # Data S1 d/dt(y182)
    d/dt(tube_gl) <- (fu_b * q_gfr * blood_gl - q_gfr * tube_gl) / v_tube_gl    # Data S1 d/dt(y183)
    d/dt(blood_pt) <- (q_renal_gl_pt * blood_gl - q_renal_pt_dt * blood_pt +
      fu_r * ps_dif_pro_ce_ve * cell_pt -
      fu_b * ps_dif_pro_ve_ce * blood_pt) / v_blood_pt                         # Data S1 d/dt(y184)
    d/dt(cell_pt) <- (-fu_r * ps_dif_pro_ce_ve * cell_pt +
      fu_b * ps_dif_pro_ve_ce * blood_pt -
      fu_r * (cl_pgp_renal * km_pgp / (km_pgp + fu_r * cell_pt) + ps_dif_pro_ce_pt) * cell_pt +
      ps_dif_pro_pt_ce * tube_pt) / v_cell_pt                                  # Data S1 d/dt(y185)
    d/dt(tube_pt) <- (q_gfr * tube_gl - q_urine_pt_dt * tube_pt +
      fu_r * (cl_pgp_renal * km_pgp / (km_pgp + fu_r * cell_pt) + ps_dif_pro_ce_pt) * cell_pt -
      ps_dif_pro_pt_ce * tube_pt) / v_tube_pt                                  # Data S1 d/dt(y186)
    d/dt(blood_dt) <- (q_renal_pt_dt * blood_pt - q_renal_dt_cd * blood_dt +
      fu_r * ps_dif_dis_ce_ve * cell_dt -
      fu_b * ps_dif_dis_ve_ce * blood_dt) / v_blood_dt                         # Data S1 d/dt(y187)
    d/dt(cell_dt) <- (-fu_r * ps_dif_dis_ce_ve * cell_dt +
      fu_b * ps_dif_dis_ve_ce * blood_dt -
      fu_r * ps_dif_dis_ce_dt * cell_dt +
      ps_dif_dis_dt_ce * tube_dt) / v_cell_dt                                  # Data S1 d/dt(y188)
    d/dt(tube_dt) <- (q_urine_pt_dt * tube_pt - q_urine_dt_cd * tube_dt +
      fu_r * ps_dif_dis_ce_dt * cell_dt -
      ps_dif_dis_dt_ce * tube_dt) / v_tube_dt                                  # Data S1 d/dt(y189)
    d/dt(blood_cd) <- (q_renal_dt_cd * blood_dt - (q_renal - q_urine_cd) * blood_cd +
      fu_r * ps_dif_col_ce_ve * cell_cd -
      fu_b * ps_dif_col_ve_ce * blood_cd) / v_blood_cd                         # Data S1 d/dt(y190)
    d/dt(cell_cd) <- (-fu_r * ps_dif_col_ce_ve * cell_cd +
      fu_b * ps_dif_col_ve_ce * blood_cd -
      fu_r * ps_dif_col_ce_cd * cell_cd +
      ps_dif_col_cd_ce * tube_cd) / v_cell_cd                                  # Data S1 d/dt(y191)
    d/dt(tube_cd) <- (q_urine_dt_cd * tube_dt - q_urine_cd * tube_cd +
      fu_r * ps_dif_col_ce_cd * cell_cd -
      ps_dif_col_cd_ce * tube_cd) / v_tube_cd                                  # Data S1 d/dt(y192)

    #    Excretion, metabolite and exposure accumulators (amounts).
    d/dt(a_urine) <- q_urine_cd * tube_cd                                      # Data S1 d/dt(y33)
    d/dt(a_feces) <- alpha_transit * k_feces_ileum * v_ileum_lumen * ileum_lumen   # Data S1 d/dt(y29)
    d/dt(a_metab_gut) <- met_duo + met_jej + met_ile                           # Data S1 d/dt(y30)
    d/dt(a_metab_liver) <- met_liver_flux                                      # Data S1 d/dt(y32)
    d/dt(auc_blood) <- central                                                 # Data S1 d/dt(y31)

    #    Cumulative canalicular secretion into bile, i.e. the flux
    #    leaving the hepatocytes. This is the state Data S1 names
    #    `BiledDuo`.
    d/dt(a_bile_canalicular) <- bile_flux                                      # Data S1 d/dt(BiledDuo)

    #    Cumulative delivery of bile INTO THE DUODENUM. This
    #    accumulator is NOT in Data S1; it is added here because it is
    #    the quantity the paper actually reports and plots (Figure 3B
    #    bile facet, "the cumulative secretion of apixaban into the
    #    duodenum compartment"; 0.554 [0.271-1.61] umol simulated
    #    versus 0.366 umol observed at 8 h after 20 mg po). It lags
    #    a_bile_canalicular by the three-compartment EHC transit and
    #    adds no dynamics of its own.
    d/dt(a_bile_duodenum) <- bile_to_duodenum

    # -----------------------------------------------------------------
    # 5. Outputs. `central` already is the TOTAL BLOOD concentration:
    #    Methods 2.2 converts every published plasma observation to
    #    blood with Rb = 1.09 before fitting, so the model's natural
    #    output scale is blood. `Cp` recovers the plasma concentration
    #    that a clinical assay would report.
    # -----------------------------------------------------------------
    Cc <- central
    Cp <- central / rb
    met_total <- a_metab_gut + a_metab_liver              # Data S1 met_total = y30 + y32
    Cc ~ prop(propSd)
  })
}
