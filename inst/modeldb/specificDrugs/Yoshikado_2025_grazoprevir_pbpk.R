Yoshikado_2025_grazoprevir_pbpk <- function() {
  description <- paste(
    "PBPK (semi-mechanistic, five-unit tandem dispersion liver,",
    "segregated-flow intestine, enterohepatic circulation). Grazoprevir",
    "after single oral doses of 50-800 mg in HCV-infected adults, fitted",
    "with the Cluster Gauss-Newton method (CGNM) using in vitro Km values",
    "as penalised reference values (the paper's 'middle-out' analysis).",
    "The nonlinear pharmacokinetics are explained by four saturable",
    "processes: hepatic uptake by OATP1B1 and OATP1B3 (Michaelis-Menten,",
    "driven by the unbound sinusoidal concentration, plus passive",
    "diffusion), hepatic CYP3A metabolism and P-gp biliary excretion",
    "(Michaelis-Menten, driven by the unbound hepatocyte concentration),",
    "and intestinal CYP3A metabolism and apical P-gp efflux in the",
    "duodenum, jejunum and ileum enterocytes. The liver is five sinusoidal",
    "sub-compartments in series, each exchanging with its own hepatocyte",
    "sub-compartment; the intestine has a lumen, enterocyte and",
    "mucosal-blood compartment per segment plus a serosal tissue",
    "compartment; muscle, skin and adipose are perfusion-limited. An oral",
    "dose enters the stomach. Biliary P-gp output reaches the duodenal",
    "lumen through a three-compartment transit chain. The SLCO1B1",
    "c.521T>C homozygous variant scales OATP1B1 activity to 34.4% of wild",
    "type. The 10 parameters without fixed() are the CGNM medians of",
    "main-text Table 2; all other values are fixed constants from Table S1",
    "and the authors' deposited model code. CGNM is a fixed-effects",
    "least-squares method, so the model has no between-subject variability",
    "and the paper reports no residual-error model -- the propSd term is a",
    "placeholder. See the vignette Errata."
  )
  reference <- paste(
    "Yoshikado T, Aoki Y, Nakamura R, Shida S, Sugiyama Y, Chiba K.",
    "Elucidating Contributions of Drug Transporters/Enzyme to Nonlinear",
    "Pharmacokinetics of Grazoprevir by PBPK Modeling With a Cluster",
    "Gauss-Newton Method.",
    "CPT Pharmacometrics Syst Pharmacol. 2025;14(4):770-780.",
    "doi:10.1002/psp4.13314.",
    "The ODE system and the fixed constants not listed in Table S1 (the",
    "per-kg organ volumes and blood flows, the gut transit scaling alpha)",
    "are transcribed from the authors' deposited rxode2 model code",
    "(Supporting Information R script, reproduced in Data S1 'Model",
    "equations'). The estimated parameters are the median column of",
    "main-text Table 2; the fixed compound and physiological constants are",
    "Table S1A.",
    sep = " "
  )
  vignette <- "Yoshikado_2025_grazoprevir_pbpk"

  # The deposited code numbers its states GZR_1..GZR_30 plus GZR_AUC; the
  # mapping to the names below is given on every d/dt() line. The stems
  # follow this laboratory's registered models of the same framework
  # (Tsuchitani_2024_telmisartan_pbpk.R, Toshimoto_2017_irinotecan_pbpk.R,
  # Aoki_2024_bosentan_pbpk.R): `is_liver<n>` is the hepatic extracellular
  # (sinusoidal) space and `int_liver<n>` the hepatocyte space of
  # dispersion segment n; `<segment>_lumen` / `<segment>_ent` /
  # `<segment>_muc` are the lumen, enterocyte and mucosal-blood
  # compartments of the segregated-flow intestine; `serosa` is the gut
  # serosal tissue; `ehc<n>` the biliary transit chain; `auc_blood` the
  # cumulative blood AUC. `a_metabolized_gut` accumulates the amount
  # metabolised by intestinal CYP3A only (hepatic CYP3A metabolism is not
  # accumulated by the deposited code), so the canonical whole-body
  # `a_metabolized` is deliberately not reused.
  paper_specific_compartments <- c(
    "is_liver1",
    "is_liver2",
    "is_liver3",
    "is_liver4",
    "is_liver5",
    "int_liver1",
    "int_liver2",
    "int_liver3",
    "int_liver4",
    "int_liver5",
    "duodenum_lumen",
    "jejunum_lumen",
    "ileum_lumen",
    "duodenum_ent",
    "jejunum_ent",
    "ileum_ent",
    "duodenum_muc",
    "jejunum_muc",
    "ileum_muc",
    "serosa",
    "ehc1",
    "ehc2",
    "ehc3",
    "a_metabolized_gut",
    "auc_blood"
  )

  # Time in hours. The deposited code works in molar units: the blood,
  # tissue, liver, lumen and enterocyte states hold concentrations in
  # nmol/L (nM, the unit of every Km in Table 2), while the stomach, bile
  # chain, faeces and gut-metabolism states hold amounts in nmol. Doses
  # are therefore supplied in nmol. The deposited event tables convert
  # 1 mg of grazoprevir to 1274 nmol (50 mg = 63,700 nmol; 187 mg =
  # 238,238 nmol), i.e. a molecular weight of 784.9 g/mol, the
  # monohydrate; the vignette uses the same factor.
  units <- list(time = "h", dosing = "nmol", concentration = "nmol/L")
  dosing <- "stomach"

  covariateData <- list(
    SNP_SLCO1B1_RS4149056_HOM = list(
      description = "SLCO1B1 c.521T>C (rs4149056) homozygous-variant (521C/C) genotype; 1 = yes, 0 = no",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (wild type, c.521 TT)",
      notes = paste(
        "Yoshikado 2025 Methods: the OATP1B1 contribution f_OATP1B1 is",
        "estimated from the observed 1.43-fold AUC increase in c.521 CC vs",
        "TT subjects (Guo 2019), assuming the variant-to-wild-type OATP1B1",
        "activity ratio RA_mt = 0.344 (Tomita 2013; Table S1A). The",
        "deposited code simulates the CC subject as a duplicate state set",
        "(MT_GZR_*) in which only the OATP1B1-mediated hepatic uptake term",
        "is multiplied by RA_mt; this covariate reproduces that exactly.",
        "Heterozygotes are not modelled by the paper; with this indicator",
        "at 0 a subject is simulated as wild type."
      ),
      source_name = "SLCO1B1 c.521 CC (MT_GZR states, RA_mt)"
    )
  )

  # What each ODE state holds, in what amount units, in what matrix.
  compartmentData <- list(
    central = list(analyte = "grazoprevir", units = "nmol/L", specimen = "whole blood", verified = TRUE),
    is_liver1 = list(analyte = "grazoprevir", units = "nmol/L", specimen = "tissue", verified = TRUE),
    is_liver2 = list(analyte = "grazoprevir", units = "nmol/L", specimen = "tissue", verified = TRUE),
    is_liver3 = list(analyte = "grazoprevir", units = "nmol/L", specimen = "tissue", verified = TRUE),
    is_liver4 = list(analyte = "grazoprevir", units = "nmol/L", specimen = "tissue", verified = TRUE),
    is_liver5 = list(analyte = "grazoprevir", units = "nmol/L", specimen = "tissue", verified = TRUE),
    int_liver1 = list(analyte = "grazoprevir", units = "nmol/L", specimen = "tissue", verified = TRUE),
    int_liver2 = list(analyte = "grazoprevir", units = "nmol/L", specimen = "tissue", verified = TRUE),
    int_liver3 = list(analyte = "grazoprevir", units = "nmol/L", specimen = "tissue", verified = TRUE),
    int_liver4 = list(analyte = "grazoprevir", units = "nmol/L", specimen = "tissue", verified = TRUE),
    int_liver5 = list(analyte = "grazoprevir", units = "nmol/L", specimen = "tissue", verified = TRUE),
    duodenum_lumen = list(analyte = "grazoprevir", units = "nmol/L", specimen = "not applicable", verified = TRUE),
    jejunum_lumen = list(analyte = "grazoprevir", units = "nmol/L", specimen = "not applicable", verified = TRUE),
    ileum_lumen = list(analyte = "grazoprevir", units = "nmol/L", specimen = "not applicable", verified = TRUE),
    muscle = list(analyte = "grazoprevir", units = "nmol/L", specimen = "tissue", verified = TRUE),
    skin = list(analyte = "grazoprevir", units = "nmol/L", specimen = "tissue", verified = TRUE),
    adipose = list(analyte = "grazoprevir", units = "nmol/L", specimen = "tissue", verified = TRUE),
    ehc1 = list(analyte = "grazoprevir", units = "nmol", specimen = "bile", verified = TRUE),
    ehc2 = list(analyte = "grazoprevir", units = "nmol", specimen = "bile", verified = TRUE),
    ehc3 = list(analyte = "grazoprevir", units = "nmol", specimen = "bile", verified = TRUE),
    duodenum_ent = list(analyte = "grazoprevir", units = "nmol/L", specimen = "tissue", verified = TRUE),
    duodenum_muc = list(analyte = "grazoprevir", units = "nmol/L", specimen = "whole blood", verified = TRUE),
    jejunum_ent = list(analyte = "grazoprevir", units = "nmol/L", specimen = "tissue", verified = TRUE),
    jejunum_muc = list(analyte = "grazoprevir", units = "nmol/L", specimen = "whole blood", verified = TRUE),
    ileum_ent = list(analyte = "grazoprevir", units = "nmol/L", specimen = "tissue", verified = TRUE),
    ileum_muc = list(analyte = "grazoprevir", units = "nmol/L", specimen = "whole blood", verified = TRUE),
    serosa = list(analyte = "grazoprevir", units = "nmol/L", specimen = "tissue", verified = TRUE),
    stomach = list(analyte = "grazoprevir", units = "nmol", specimen = "not applicable", verified = TRUE),
    a_feces = list(analyte = "grazoprevir", units = "nmol", specimen = "faeces", verified = TRUE),
    a_metabolized_gut = list(analyte = "grazoprevir", units = "nmol", specimen = "not applicable", verified = TRUE),
    auc_blood = list(analyte = "grazoprevir", units = "nmol/L*h", specimen = "whole blood", verified = TRUE)
  )

  population <- list(
    species = "human",
    dose_range = "single oral doses of 50, 100, 200, 400, 600 and 800 mg",
    disease_state = "chronic hepatitis C (genotype 1 and genotype 3 HCV-infected patients)",
    notes = paste(
      "Yoshikado 2025 fits the mean grazoprevir concentration-time",
      "profiles of a dose-escalation study (50-800 mg) reported only as a",
      "conference abstract (Caro et al., AASLD 2011, MK-5172 in genotype 1",
      "and 3 HCV-infected patients), together with the cumulative faecal",
      "excretion of 187 mg from the Zepatier label and the 1.43-fold AUC",
      "ratio in SLCO1B1 c.521 CC vs TT subjects (Guo 2019). Subject counts",
      "and demographics are not reported. Body weight is fixed at 70 kg",
      "because the clinical source does not report it (Table S1A); every",
      "organ volume and blood flow in the model is that per-kg value times",
      "70 kg."
    )
  )

  ini({
    # =================================================================
    # ESTIMATED PARAMETERS -- main-text Table 2, 'Median' column: the
    # PBPK-CGNM analysis with in vitro Km values as penalised references
    # and optimised weights. With the authors' deposited code and data,
    # this median vector gives SSR' = 2.537, matching the minimum SSR'
    # of 2.54 printed in Figure 3A, so the medians reproduce the
    # accepted fit. Km values are unbound concentrations in nM.
    # =================================================================
    beta_hep <- 0.658
    label("Fraction of hepatocyte disposition that is elimination rather than sinusoidal efflux (unitless)") # Table 2 beta, median 0.658
    lcl_int_all <- log(12900)
    label("Overall hepatic intrinsic clearance, unbound (L/h)") # Table 2 CLint,all, median 12,900
    f_bile <- 0.882
    label("Fraction of hepatocyte elimination that is P-gp biliary excretion rather than CYP3A metabolism (unitless)") # Table 2 fbile, median 0.882
    f_oatp1b1 <- 0.197
    label("Fraction of hepatic active uptake mediated by OATP1B1 (unitless)") # Table 2 fOATP1B1, median 0.197
    lkm_oatp1b1 <- log(136)
    label("Michaelis-Menten constant of OATP1B1 uptake, unbound (nmol/L)") # Table 2 Km,OATP1B1, median 136 nM
    lkm_cyp3a <- log(1570)
    label("Michaelis-Menten constant of CYP3A metabolism, unbound (nmol/L)") # Table 2 Km,CYP3A, median 1570 nM
    lkm_pgp <- log(347)
    label("Michaelis-Menten constant of P-gp efflux, unbound (nmol/L)") # Table 2 Km,P-gp, median 347 nM
    lumen_ratio <- 0.668
    label("Lumen-to-enterocyte relative to basolateral passive diffusion clearance, PSdif,ent,LE/PSdif,ent,BE (unitless)") # Table 2 LR, median 0.668
    lps_dif_ent_be <- log(3.86)
    label("Passive diffusion clearance across the enterocyte basolateral membrane, blood to enterocyte (L/h)") # Table 2 PSdif,ent,BE, median 3.86
    lkpuu <- log(71.4)
    label("Liver-to-blood unbound concentration ratio Kp,uu (unitless)") # Table 2 Kp,uu, median 71.4

    # =================================================================
    # FIXED COMPOUND CONSTANTS -- Table S1A (and the identical
    # assignments in the deposited code).
    # =================================================================
    r_km <- fixed(0.235)
    label("Ratio of the OATP1B3 to the OATP1B1 Michaelis-Menten constant (unitless)") # code rKm=0.235; in vitro 39.9/170 nM (Results); Table 2 footnote prints 0.236
    gamma_hep <- fixed(0.24)
    label("Hepatic passive influx-to-efflux clearance ratio PSdif,inf/PSeff (unitless)") # Table S1A gamma 0.24 (in silico prediction)
    gamma_ent <- fixed(0.24)
    label("Enterocyte basolateral-to-apical passive diffusion clearance ratio (unitless)") # Table S1A gamma,ent 0.24, assumed equal to liver
    apical_ratio <- fixed(1)
    label("Enterocyte-to-lumen relative to enterocyte-to-blood passive diffusion clearance, PSdif,ent,EL/PSdif,ent,EB (unitless)") # Table S1A AR 1
    fu_b <- fixed(0.0241)
    label("Unbound fraction of grazoprevir in blood (unitless)") # Table S1A fB 0.0241
    fu_h <- fixed(0.0173)
    label("Unbound fraction of grazoprevir in hepatocytes (unitless)") # Table S1A fh 0.0173, measured in plated hepatocytes
    fu_gut <- fixed(0.0197)
    label("Unbound fraction of grazoprevir in enterocytes (unitless)") # Table S1A fgut 0.0197
    lkp_adipose <- fixed(log(0.126))
    label("Adipose-to-blood partition coefficient (unitless)") # Table S1A Kp,adipose 0.126
    lkp_gut <- fixed(log(0.26))
    label("Gut-serosa-to-blood partition coefficient (unitless)") # Table S1A Kp,gut 0.26
    lkp_muscle <- fixed(log(0.107))
    label("Muscle-to-blood partition coefficient (unitless)") # Table S1A Kp,muscle 0.107
    lkp_skin <- fixed(log(0.446))
    label("Skin-to-blood partition coefficient (unitless)") # Table S1A Kp,skin 0.446
    sf_kp <- fixed(1)
    label("Scaling factor applied to every tissue-to-blood partition coefficient (unitless)") # Table S1A SF_Kp 1
    cl_r <- fixed(0)
    label("Renal clearance of grazoprevir from blood (L/h)") # Table S1A CLr 0
    ra_oatp1b1_mt <- fixed(0.344)
    label("OATP1B1 activity of the c.521C variant relative to wild type (unitless)") # Table S1A RA OATP1B1_mt/wt 0.344
    vmax_ratio_cyp3a <- fixed(0.00741)
    label("Intestinal-to-hepatic CYP3A activity ratio (unitless)") # Table S1A RA CYP3A_gut/liver 0.00741
    vmax_ratio_pgp <- fixed(0.0365)
    label("Intestinal-to-hepatic P-gp activity ratio (unitless)") # Table S1A RA P-gp_gut/liver 0.0365
    k_bile <- fixed(0.605)
    label("Biliary transit rate constant (1/h)") # Table S1A kbile 0.605
    k_stomach <- fixed(2.5)
    label("Gastric emptying rate constant (1/h)") # Table S1A kstomach 2.5
    kfeces_duodenum <- fixed(5.12)
    label("Duodenal lumen transit rate constant (1/h)") # Table S1A kfeces_duodenum 5.12
    kfeces_jejunum <- fixed(1.47)
    label("Jejunal lumen transit rate constant (1/h)") # Table S1A kfeces_jejunum 1.47
    kfeces_ileum <- fixed(1.21)
    label("Ileal lumen transit rate constant to faeces (1/h)") # Table S1A kfeces_ileum 1.21
    alpha_transit <- fixed(1)
    label("Scaling factor applied to every intestinal transit rate constant (unitless)") # deposited code Alpha=1 (not listed in Table S1)

    # Regional distribution of intestinal processes over duodenum,
    # jejunum and ileum. The P-gp fractions are the deposited code's
    # (0.057 / 0.273 / 0.670); Table S1A prints 0.0197 / 0.371 / 0.61.
    # The two sets give SSR' 2.5367 and 2.5364 at the Table 2 medians,
    # so the choice is immaterial to the fit; the code is used because it
    # is the model that produced Table 2.
    f_cyp_duodenum <- fixed(0.138)
    label("Fraction of intestinal CYP3A in the duodenum (unitless)") # Table S1A Fcyp_duodenum 0.138
    f_cyp_jejunum <- fixed(0.545)
    label("Fraction of intestinal CYP3A in the jejunum (unitless)") # Table S1A Fcyp_jejunum 0.545
    f_cyp_ileum <- fixed(0.317)
    label("Fraction of intestinal CYP3A in the ileum (unitless)") # Table S1A Fcyp_ileum 0.317
    f_dif_duodenum <- fixed(0.076)
    label("Fraction of intestinal passive diffusion clearance in the duodenum (unitless)") # Table S1A Fdif_duodenum 0.076
    f_dif_jejunum <- fixed(0.395)
    label("Fraction of intestinal passive diffusion clearance in the jejunum (unitless)") # Table S1A Fdif_jejunum 0.395
    f_dif_ileum <- fixed(0.529)
    label("Fraction of intestinal passive diffusion clearance in the ileum (unitless)") # Table S1A Fdif_ileum 0.529
    f_pgp_duodenum <- fixed(0.057)
    label("Fraction of intestinal P-gp in the duodenum (unitless)") # deposited code Fpgp_duodenum=0.057 (Table S1A prints 0.0197)
    f_pgp_jejunum <- fixed(0.273)
    label("Fraction of intestinal P-gp in the jejunum (unitless)") # deposited code Fpgp_jejunum=0.273 (Table S1A prints 0.371)
    f_pgp_ileum <- fixed(0.670)
    label("Fraction of intestinal P-gp in the ileum (unitless)") # deposited code Fpgp_ileum=0.670 (Table S1A prints 0.61)

    # =================================================================
    # FIXED PHYSIOLOGY -- per kg of body weight. Table S1A states that
    # these follow the standard protocol of Yoshikado 2016 (its ref. 2)
    # without listing them; the values are the coefficients of BW in the
    # deposited code.
    # =================================================================
    bw <- fixed(70)
    label("Body weight (kg)") # Table S1A BW 70, assumed
    v_central_kg <- fixed(0.0743)
    label("Central (blood) volume per kg (L/kg)") # deposited code 0.0743*BW
    v_is_liver_kg <- fixed(0.0067)
    label("Hepatic extracellular volume per kg (L/kg)") # deposited code 0.0067*BW, split over 5 segments
    v_int_liver_kg <- fixed(0.0174)
    label("Hepatocyte volume per kg (L/kg)") # deposited code 0.0174*BW, split over 5 segments
    v_muscle_kg <- fixed(0.429)
    label("Muscle volume per kg (L/kg)") # deposited code 0.429*BW
    v_skin_kg <- fixed(0.111)
    label("Skin volume per kg (L/kg)") # deposited code 0.111*BW
    v_adipose_kg <- fixed(0.143)
    label("Adipose volume per kg (L/kg)") # deposited code 0.143*BW
    v_serosa_kg <- fixed(0.00702)
    label("Gut serosa volume per kg (L/kg)") # deposited code 0.00702*BW
    v_duodenum_lumen_kg <- fixed(0.00093)
    label("Duodenal lumen volume per kg (L/kg)") # deposited code 0.00093*BW
    v_jejunum_lumen_kg <- fixed(0.00331)
    label("Jejunal lumen volume per kg (L/kg)") # deposited code 0.00331*BW
    v_ileum_lumen_kg <- fixed(0.00395)
    label("Ileal lumen volume per kg (L/kg)") # deposited code 0.00395*BW
    v_duodenum_ent_kg <- fixed(0.00084)
    label("Duodenal enterocyte volume per kg (L/kg)") # deposited code 0.00084*BW
    v_jejunum_ent_kg <- fixed(0.00299)
    label("Jejunal enterocyte volume per kg (L/kg)") # deposited code 0.00299*BW
    v_ileum_ent_kg <- fixed(0.00357)
    label("Ileal enterocyte volume per kg (L/kg)") # deposited code 0.00357*BW
    v_duodenum_muc_kg <- fixed(0.000113)
    label("Duodenal mucosal blood volume per kg (L/kg)") # deposited code 0.000113*BW
    v_jejunum_muc_kg <- fixed(0.000400)
    label("Jejunal mucosal blood volume per kg (L/kg)") # deposited code 0.000400*BW
    v_ileum_muc_kg <- fixed(0.000477)
    label("Ileal mucosal blood volume per kg (L/kg)") # deposited code 0.000477*BW
    q_hepatic_kg <- fixed(1.24)
    label("Hepatic blood flow per kg (L/h/kg)") # deposited code 1.24*BW
    q_muscle_kg <- fixed(0.642)
    label("Muscle blood flow per kg (L/h/kg)") # deposited code 0.642*BW
    q_skin_kg <- fixed(0.257)
    label("Skin blood flow per kg (L/h/kg)") # deposited code 0.257*BW
    q_adipose_kg <- fixed(0.223)
    label("Adipose blood flow per kg (L/h/kg)") # deposited code 0.223*BW
    q_serosa_kg <- fixed(0.744)
    label("Gut serosal blood flow per kg (L/h/kg)") # deposited code 0.744*BW
    q_duodenum_muc_kg <- fixed(0.0164)
    label("Duodenal mucosal blood flow per kg (L/h/kg)") # deposited code 0.0164*BW
    q_jejunum_muc_kg <- fixed(0.0900)
    label("Jejunal mucosal blood flow per kg (L/h/kg)") # deposited code 0.0900*BW
    q_ileum_muc_kg <- fixed(0.0796)
    label("Ileal mucosal blood flow per kg (L/h/kg)") # deposited code 0.0796*BW

    # =================================================================
    # RESIDUAL ERROR. CGNM minimises an unweighted sum of squared log10
    # residuals plus penalty terms (main-text Equation 2) and estimates
    # no residual-error model. The term below exists only so the model is
    # a valid nlmixr2 object; it is not an estimate from the paper. Same
    # convention as Tsuchitani_2024_telmisartan_pbpk.R.
    # =================================================================
    propSd <- fixed(0.10)
    label("Proportional residual error placeholder (fraction)") # not reported by Yoshikado 2025; placeholder only
  })

  model({
    # -----------------------------------------------------------------
    # 1. Back-transform and scale to the 70 kg body.
    # -----------------------------------------------------------------
    cl_int_all <- exp(lcl_int_all)
    km_oatp1b1 <- exp(lkm_oatp1b1)
    km_cyp3a <- exp(lkm_cyp3a)
    km_pgp <- exp(lkm_pgp)
    ps_dif_ent_be <- exp(lps_dif_ent_be)
    kpuu <- exp(lkpuu)
    kp_adipose <- exp(lkp_adipose)
    kp_gut <- exp(lkp_gut)
    kp_muscle <- exp(lkp_muscle)
    kp_skin <- exp(lkp_skin)

    v_central <- v_central_kg * bw
    v_is_seg <- 0.2 * v_is_liver_kg * bw
    v_int_seg <- 0.2 * v_int_liver_kg * bw
    v_muscle <- v_muscle_kg * bw
    v_skin <- v_skin_kg * bw
    v_adipose <- v_adipose_kg * bw
    v_serosa <- v_serosa_kg * bw
    v_duodenum_lumen <- v_duodenum_lumen_kg * bw
    v_jejunum_lumen <- v_jejunum_lumen_kg * bw
    v_ileum_lumen <- v_ileum_lumen_kg * bw
    v_duodenum_ent <- v_duodenum_ent_kg * bw
    v_jejunum_ent <- v_jejunum_ent_kg * bw
    v_ileum_ent <- v_ileum_ent_kg * bw
    v_duodenum_muc <- v_duodenum_muc_kg * bw
    v_jejunum_muc <- v_jejunum_muc_kg * bw
    v_ileum_muc <- v_ileum_muc_kg * bw
    q_hepatic <- q_hepatic_kg * bw
    q_muscle <- q_muscle_kg * bw
    q_skin <- q_skin_kg * bw
    q_adipose <- q_adipose_kg * bw
    q_serosa <- q_serosa_kg * bw
    q_duodenum_muc <- q_duodenum_muc_kg * bw
    q_jejunum_muc <- q_jejunum_muc_kg * bw
    q_ileum_muc <- q_ileum_muc_kg * bw
    # Hepatic artery: hepatic outflow minus the portal (serosal and
    # mucosal) inflows.
    q_ha <- q_hepatic - q_serosa - (q_duodenum_muc + q_jejunum_muc + q_ileum_muc)

    # -----------------------------------------------------------------
    # 2. Hybrid-to-elementary hepatic parameters (Yoshikado 2016 /
    #    Toshimoto 2017 framework, as written in the deposited code).
    #    Rdif = PSdif,inf / PSact,inf follows from Kp,uu, gamma and beta;
    #    CLint,all = (PSact,inf + PSdif,inf) * beta.
    # -----------------------------------------------------------------
    r_dif <- gamma_hep * (1 - beta_hep) / (kpuu - gamma_hep * (1 - beta_hep))
    ps_act_inf <- cl_int_all / beta_hep / (1 + r_dif)
    ps_dif_inf <- r_dif / (1 + r_dif) * cl_int_all / beta_hep
    ps_dif_eff <- ps_dif_inf / gamma_hep
    cl_int_elim <- cl_int_all / (1 - beta_hep) * r_dif / (1 + r_dif) / gamma_hep

    # OATP1B1 activity is reduced in c.521 CC subjects (the deposited
    # code's MT_GZR states multiply only this term by RA_mt).
    ra_oatp1b1 <- 1 - SNP_SLCO1B1_RS4149056_HOM * (1 - ra_oatp1b1_mt)
    km_oatp1b3 <- r_km * km_oatp1b1
    vmax_oatp1b1 <- ra_oatp1b1 * km_oatp1b1 * f_oatp1b1 * ps_act_inf
    vmax_oatp1b3 <- km_oatp1b3 * (1 - f_oatp1b1) * ps_act_inf
    vmax_pgp <- km_pgp * f_bile * cl_int_elim
    vmax_cyp3a <- km_cyp3a * (1 - f_bile) * cl_int_elim
    vmax_pgp_ent <- vmax_pgp * vmax_ratio_pgp
    vmax_cyp3a_ent <- vmax_cyp3a * vmax_ratio_cyp3a
    ps_dif_ent_eb <- ps_dif_ent_be / gamma_ent

    # Blood concentrations leaving the perfusion-limited tissues.
    cb_muscle <- muscle / kp_muscle / sf_kp
    cb_skin <- skin / kp_skin / sf_kp
    cb_adipose <- adipose / kp_adipose / sf_kp
    cb_serosa <- serosa / kp_gut / sf_kp

    # -----------------------------------------------------------------
    # 3. Hepatic transport and elimination per dispersion segment (each
    #    segment carries one fifth of every clearance). Uptake is
    #    saturable OATP1B1 + OATP1B3 plus passive diffusion on the unbound
    #    sinusoidal concentration; elimination is saturable P-gp biliary
    #    excretion + CYP3A metabolism on the unbound hepatocyte
    #    concentration.
    # -----------------------------------------------------------------
    uptake1 <- 0.2 * fu_b * (vmax_oatp1b1 / (km_oatp1b1 + fu_b * is_liver1) + vmax_oatp1b3 / (km_oatp1b3 + fu_b * is_liver1) + ps_dif_inf) * is_liver1
    uptake2 <- 0.2 * fu_b * (vmax_oatp1b1 / (km_oatp1b1 + fu_b * is_liver2) + vmax_oatp1b3 / (km_oatp1b3 + fu_b * is_liver2) + ps_dif_inf) * is_liver2
    uptake3 <- 0.2 * fu_b * (vmax_oatp1b1 / (km_oatp1b1 + fu_b * is_liver3) + vmax_oatp1b3 / (km_oatp1b3 + fu_b * is_liver3) + ps_dif_inf) * is_liver3
    uptake4 <- 0.2 * fu_b * (vmax_oatp1b1 / (km_oatp1b1 + fu_b * is_liver4) + vmax_oatp1b3 / (km_oatp1b3 + fu_b * is_liver4) + ps_dif_inf) * is_liver4
    uptake5 <- 0.2 * fu_b * (vmax_oatp1b1 / (km_oatp1b1 + fu_b * is_liver5) + vmax_oatp1b3 / (km_oatp1b3 + fu_b * is_liver5) + ps_dif_inf) * is_liver5

    efflux1 <- 0.2 * fu_h * ps_dif_eff * int_liver1
    efflux2 <- 0.2 * fu_h * ps_dif_eff * int_liver2
    efflux3 <- 0.2 * fu_h * ps_dif_eff * int_liver3
    efflux4 <- 0.2 * fu_h * ps_dif_eff * int_liver4
    efflux5 <- 0.2 * fu_h * ps_dif_eff * int_liver5

    elim1 <- 0.2 * fu_h * (vmax_pgp / (km_pgp + fu_h * int_liver1) + vmax_cyp3a / (km_cyp3a + fu_h * int_liver1)) * int_liver1
    elim2 <- 0.2 * fu_h * (vmax_pgp / (km_pgp + fu_h * int_liver2) + vmax_cyp3a / (km_cyp3a + fu_h * int_liver2)) * int_liver2
    elim3 <- 0.2 * fu_h * (vmax_pgp / (km_pgp + fu_h * int_liver3) + vmax_cyp3a / (km_cyp3a + fu_h * int_liver3)) * int_liver3
    elim4 <- 0.2 * fu_h * (vmax_pgp / (km_pgp + fu_h * int_liver4) + vmax_cyp3a / (km_cyp3a + fu_h * int_liver4)) * int_liver4
    elim5 <- 0.2 * fu_h * (vmax_pgp / (km_pgp + fu_h * int_liver5) + vmax_cyp3a / (km_cyp3a + fu_h * int_liver5)) * int_liver5

    # Hepatic P-gp biliary excretion (feeds the bile chain).
    pgp_bile <- 0.2 * fu_h * vmax_pgp *
      (int_liver1 / (km_pgp + fu_h * int_liver1) +
        int_liver2 / (km_pgp + fu_h * int_liver2) +
        int_liver3 / (km_pgp + fu_h * int_liver3) +
        int_liver4 / (km_pgp + fu_h * int_liver4) +
        int_liver5 / (km_pgp + fu_h * int_liver5))

    # -----------------------------------------------------------------
    # 4. Enterocyte fluxes per segment (nmol/h). Apical secretion back to
    #    the lumen is saturable P-gp plus passive diffusion; gut CYP3A
    #    metabolism is saturable.
    # -----------------------------------------------------------------
    secr_duodenum <- fu_gut * (vmax_pgp_ent * f_pgp_duodenum / (km_pgp + fu_gut * duodenum_ent) + apical_ratio * ps_dif_ent_eb * f_dif_duodenum) * duodenum_ent
    secr_jejunum <- fu_gut * (vmax_pgp_ent * f_pgp_jejunum / (km_pgp + fu_gut * jejunum_ent) + apical_ratio * ps_dif_ent_eb * f_dif_jejunum) * jejunum_ent
    secr_ileum <- fu_gut * (vmax_pgp_ent * f_pgp_ileum / (km_pgp + fu_gut * ileum_ent) + apical_ratio * ps_dif_ent_eb * f_dif_ileum) * ileum_ent

    met_duodenum <- fu_gut * vmax_cyp3a_ent * f_cyp_duodenum / (km_cyp3a + fu_gut * duodenum_ent) * duodenum_ent
    met_jejunum <- fu_gut * vmax_cyp3a_ent * f_cyp_jejunum / (km_cyp3a + fu_gut * jejunum_ent) * jejunum_ent
    met_ileum <- fu_gut * vmax_cyp3a_ent * f_cyp_ileum / (km_cyp3a + fu_gut * ileum_ent) * ileum_ent

    # -----------------------------------------------------------------
    # 5. ODEs -- deposited code; the code's state name ends every line.
    # -----------------------------------------------------------------
    d/dt(central) <- 1 / v_central *
      (q_hepatic * (is_liver5 - central) - cl_r * central -
        q_muscle * (central - cb_muscle) -
        q_skin * (central - cb_skin) -
        q_adipose * (central - cb_adipose)) # GZR_1

    # Liver sinusoidal segments; segment 1 receives the hepatic artery,
    # the serosal and the three mucosal-blood drainages.
    d/dt(is_liver1) <- 1 / v_is_seg *
      (q_ha * central + q_serosa * cb_serosa +
        (q_duodenum_muc * duodenum_muc + q_jejunum_muc * jejunum_muc + q_ileum_muc * ileum_muc) -
        q_hepatic * is_liver1 - uptake1 + efflux1) # GZR_2
    d/dt(is_liver2) <- 1 / v_is_seg * (q_hepatic * (is_liver1 - is_liver2) - uptake2 + efflux2) # GZR_4
    d/dt(is_liver3) <- 1 / v_is_seg * (q_hepatic * (is_liver2 - is_liver3) - uptake3 + efflux3) # GZR_6
    d/dt(is_liver4) <- 1 / v_is_seg * (q_hepatic * (is_liver3 - is_liver4) - uptake4 + efflux4) # GZR_8
    d/dt(is_liver5) <- 1 / v_is_seg * (q_hepatic * (is_liver4 - is_liver5) - uptake5 + efflux5) # GZR_10

    # Hepatocyte segments.
    d/dt(int_liver1) <- 1 / v_int_seg * (uptake1 - efflux1 - elim1) # GZR_3
    d/dt(int_liver2) <- 1 / v_int_seg * (uptake2 - efflux2 - elim2) # GZR_5
    d/dt(int_liver3) <- 1 / v_int_seg * (uptake3 - efflux3 - elim3) # GZR_7
    d/dt(int_liver4) <- 1 / v_int_seg * (uptake4 - efflux4 - elim4) # GZR_9
    d/dt(int_liver5) <- 1 / v_int_seg * (uptake5 - efflux5 - elim5) # GZR_11

    # Small-intestinal lumen; the duodenum receives gastric emptying and
    # bile. The ileum empties directly to faeces.
    d/dt(duodenum_lumen) <- 1 / v_duodenum_lumen *
      (k_stomach * stomach + k_bile * ehc3 -
        (alpha_transit * kfeces_duodenum * v_duodenum_lumen + lumen_ratio * ps_dif_ent_be * f_dif_duodenum) * duodenum_lumen +
        secr_duodenum) # GZR_12
    d/dt(jejunum_lumen) <- 1 / v_jejunum_lumen *
      (alpha_transit * kfeces_duodenum * v_duodenum_lumen * duodenum_lumen -
        (alpha_transit * kfeces_jejunum * v_jejunum_lumen + lumen_ratio * ps_dif_ent_be * f_dif_jejunum) * jejunum_lumen +
        secr_jejunum) # GZR_13
    d/dt(ileum_lumen) <- 1 / v_ileum_lumen *
      (alpha_transit * kfeces_jejunum * v_jejunum_lumen * jejunum_lumen -
        (alpha_transit * kfeces_ileum * v_ileum_lumen + lumen_ratio * ps_dif_ent_be * f_dif_ileum) * ileum_lumen +
        secr_ileum) # GZR_14

    # Perfusion-limited tissues.
    d/dt(muscle) <- 1 / v_muscle * (q_muscle * (central - cb_muscle)) # GZR_15
    d/dt(skin) <- 1 / v_skin * (q_skin * (central - cb_skin)) # GZR_16
    d/dt(adipose) <- 1 / v_adipose * (q_adipose * (central - cb_adipose)) # GZR_17

    # Bile transit chain (amounts).
    d/dt(ehc1) <- pgp_bile - k_bile * ehc1 # GZR_18
    d/dt(ehc2) <- k_bile * (ehc1 - ehc2) # GZR_19
    d/dt(ehc3) <- k_bile * (ehc2 - ehc3) # GZR_20

    # Enterocytes and mucosal blood, per segment.
    d/dt(duodenum_ent) <- 1 / v_duodenum_ent *
      (lumen_ratio * ps_dif_ent_be * f_dif_duodenum * duodenum_lumen - secr_duodenum -
        fu_gut * ps_dif_ent_eb * f_dif_duodenum * duodenum_ent +
        fu_b * ps_dif_ent_be * f_dif_duodenum * duodenum_muc - met_duodenum) # GZR_21
    d/dt(duodenum_muc) <- 1 / v_duodenum_muc *
      (q_duodenum_muc * (central - duodenum_muc) +
        fu_gut * ps_dif_ent_eb * f_dif_duodenum * duodenum_ent -
        fu_b * ps_dif_ent_be * f_dif_duodenum * duodenum_muc) # GZR_22
    d/dt(jejunum_ent) <- 1 / v_jejunum_ent *
      (lumen_ratio * ps_dif_ent_be * f_dif_jejunum * jejunum_lumen - secr_jejunum -
        fu_gut * ps_dif_ent_eb * f_dif_jejunum * jejunum_ent +
        fu_b * ps_dif_ent_be * f_dif_jejunum * jejunum_muc - met_jejunum) # GZR_23
    d/dt(jejunum_muc) <- 1 / v_jejunum_muc *
      (q_jejunum_muc * (central - jejunum_muc) +
        fu_gut * ps_dif_ent_eb * f_dif_jejunum * jejunum_ent -
        fu_b * ps_dif_ent_be * f_dif_jejunum * jejunum_muc) # GZR_24
    d/dt(ileum_ent) <- 1 / v_ileum_ent *
      (lumen_ratio * ps_dif_ent_be * f_dif_ileum * ileum_lumen - secr_ileum -
        fu_gut * ps_dif_ent_eb * f_dif_ileum * ileum_ent +
        fu_b * ps_dif_ent_be * f_dif_ileum * ileum_muc - met_ileum) # GZR_25
    d/dt(ileum_muc) <- 1 / v_ileum_muc *
      (q_ileum_muc * (central - ileum_muc) +
        fu_gut * ps_dif_ent_eb * f_dif_ileum * ileum_ent -
        fu_b * ps_dif_ent_be * f_dif_ileum * ileum_muc) # GZR_26

    # Gut serosa (non-absorptive limb of the segregated-flow intestine).
    d/dt(serosa) <- 1 / v_serosa * (q_serosa * (central - cb_serosa)) # GZR_27

    # Stomach (amount); the oral dose enters here.
    d/dt(stomach) <- -k_stomach * stomach # GZR_28

    # Cumulative faecal excretion, gut CYP3A metabolism and blood AUC.
    d/dt(a_feces) <- alpha_transit * kfeces_ileum * v_ileum_lumen * ileum_lumen # GZR_29
    d/dt(a_metabolized_gut) <- met_duodenum + met_jejunum + met_ileum # GZR_30
    d/dt(auc_blood) <- central # GZR_AUC

    # -----------------------------------------------------------------
    # 6. Observation. The model's central compartment and the paper's
    #    figures are whole-blood concentrations in nmol/L.
    # -----------------------------------------------------------------
    Cc <- central
    Cc ~ prop(propSd)
  })
}
