Clements_2026_belantamab <- function() {
  description <- paste(
    "Simultaneous population PKPD model of the antibody-drug conjugate",
    "belantamab mafodotin, soluble BCMA (sBCMA) and serum M-protein in",
    "patients with relapsed/refractory multiple myeloma (Clements 2026;",
    "DREAMM-1, -2, -3 and -5 monotherapy). Belantamab mafodotin follows a",
    "two-compartment model with linear clearance and quasi-equilibrium",
    "target-mediated drug disposition (TMDD) by sBCMA, with the",
    "drug-sBCMA complex internalised at the linear elimination rate (so",
    "total-drug clearance is linear) and total drug distributing to the",
    "peripheral compartment. sBCMA is synthesised into a reservoir in",
    "proportion to the current-to-baseline M-protein ratio and diffuses",
    "into the central total-sBCMA pool. M-protein follows logistic growth",
    "(carrying capacity 150 g/L) with a kill rate proportional to the free",
    "belantamab mafodotin concentration. Covariates: baseline albumin,",
    "sBCMA, body weight and sex on clearance; body weight and sex on",
    "central volume; baseline sBCMA and extramedullary disease on the",
    "M-protein kill rate. Outputs: total belantamab mafodotin (Cc), free",
    "and complexed sBCMA (with the paper's assay-dilution corrections), and",
    "M-protein, all in the paper's log-transformed residual error model.",
    sep = " "
  )
  reference <- paste(
    "Clements JD, Sorribes IC, Kaullen J, Neumar C, Taylor A, Yue MS,",
    "Chen X, Struemper H, Ferron-Brady G. Simultaneous Population PKPD",
    "Analysis of Belantamab Mafodotin, Soluble BCMA, and Serum M-Protein in",
    "Relapsed/Refractory Multiple Myeloma. CPT Pharmacometrics Syst",
    "Pharmacol. 2026;15(10):e70337. doi:10.1002/psp4.70337.",
    "Builds on the population PK model of Rathi 2021",
    "(modellib('Rathi_2021_belantamab')) and the M-protein model of Collins",
    "2023 (modellib('Collins_2023_belantamab_mprotein')).",
    sep = " "
  )
  vignette <- "Clements_2026_belantamab"
  # Dose in mg; the ADC state is carried in mg and converted to nM for the
  # binding equations. Every plasma output (total ADC, free and complexed
  # sBCMA) is in nM, as in the paper's dataset; M-protein is in g/L.
  units <- list(time = "day", dosing = "mg", concentration = "nM")

  covariateData <- list(
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-fixed baseline. Power model normalised to the analysis-dataset",
        "median 74.65 kg (control stream 'WTBL/74.65'; Table 1 overall",
        "median 74.7 kg), exponent 0.488 on CL and 0.506 on Vc.",
        sep = " "
      ),
      source_name = "WTBL"
    ),
    ALB = list(
      description = "Baseline serum albumin",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-fixed baseline. Power model normalised to 38 g/L (control",
        "stream 'ALBBL/38'; Table 1 overall median 38.0 g/L), exponent",
        "-0.879 on CL.",
        sep = " "
      ),
      source_name = "ALBBL"
    ),
    SBCMA = list(
      description = "Baseline serum soluble B-cell maturation antigen (sBCMA) concentration",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-fixed baseline in ng/mL. Three roles: (1) power covariate",
        "normalised to 72.85 ng/mL (control stream 'SBCMABL/72.85'; Table 1",
        "overall median 72.9 ng/mL) with exponent 0.258 on CL and 0.318 on",
        "the M-protein kill rate; (2) converted to nM with the control-stream",
        "molecular weight MWSBCMA = 6 kDa (SBCMA / 6) to give the initial",
        "total sBCMA and reservoir concentrations; (3) sets the sBCMA",
        "synthesis rate KSYN = KDEG * SBCMA / 6.",
        sep = " "
      ),
      source_name = "SBCMABL"
    ),
    SEXF = list(
      description = "Biological sex, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = paste(
        "Source column SEXN codes 1 = male and 2 = female; the control stream",
        "applies (1 + theta) when SEXN = 2, so SEXF = 1 for SEXN = 2.",
        "Proportional shifts -0.126 on CL and -0.136 on Vc for females.",
        sep = " "
      ),
      source_name = "SEXN"
    ),
    DIS_EMD = list(
      description = "Extramedullary disease at screening",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no extramedullary disease)",
      notes = paste(
        "Source column MEDFL (1 = present). The control stream maps the",
        "missing code -999 to the reference category. Proportional shift",
        "-0.684 on the M-protein kill rate (KDEATH) when DIS_EMD = 1.",
        sep = " "
      ),
      source_name = "MEDFL"
    ),
    MCPROT = list(
      description = "Baseline serum monoclonal (M) protein concentration",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-fixed baseline in g/L (register default is g/dL; 1 g/dL =",
        "10 g/L). Seeds the initial condition of the `mprotein` state",
        "(control stream A_0(4) = MPROTBL) and normalises the M-protein",
        "drive on sBCMA synthesis, KSYN * (MPROT / MPROTBL). Analysis",
        "required a quantifiable baseline above 0.5 g/L; Table 1 overall",
        "median 18.0 g/L (range 0.7-102).",
        sep = " "
      ),
      source_name = "MPROTBL"
    )
  )

  covariatesDataExcluded <- list(
    B2M = list(
      description = "Baseline serum beta-2-microglobulin",
      units = "nmol/L",
      type = "continuous",
      notes = paste(
        "Significant on CL in forward addition, then removed from the final",
        "model 'due to a small effect and large RSE' (Results). The control",
        "stream keeps the term as (B2MBL/320.34)^THETA(11) with THETA(11)",
        "fixed to 0, so it contributes a factor of exactly 1.",
        sep = " "
      )
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Tested on CL, Vc and Vp (Table S2); no significant effect (Results)."
    ),
    IGG = list(
      description = "Baseline serum immunoglobulin G",
      units = "g/L",
      type = "continuous",
      notes = paste(
        "Tested on CL, Vc, KDEATH and KGR (Table S2) and again on CL after",
        "stepwise covariate modelling (Methods); not retained.",
        sep = " "
      )
    )
  )

  compartmentData <- list(
    central = list(
      analyte = "belantamab mafodotin (total: free + sBCMA-bound)",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(analyte = "belantamab mafodotin (total)", units = "mg", specimen = "plasma", verified = TRUE),
    total_target = list(
      analyte = "soluble BCMA (total: free + belantamab mafodotin-bound)",
      units = "nM",
      specimen = "serum",
      verified = TRUE
    ),
    mprotein = list(analyte = "M-protein", units = "g/L", specimen = "serum", verified = TRUE),
    precursor1 = list(
      analyte = "soluble BCMA (synthesis reservoir)",
      units = "nM",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 510L,
    n_studies = 4L,
    age_range = "34-89 years",
    age_median = "66.0 years",
    weight_range = "37.0-142 kg",
    weight_median = "74.7 kg",
    sex_female_pct = 45.9,
    race_ethnicity = c(White = 80.4, Black = 8.0, Asian = 11.4, Missing = 0.2),
    disease_state = paste(
      "Relapsed/refractory multiple myeloma with at least one",
      "pharmacodynamic measurement (serum M-protein or sBCMA) and a",
      "quantifiable baseline serum M-protein above 0.5 g/L.",
      sep = " "
    ),
    dose_range = paste(
      "Belantamab mafodotin monotherapy by IV infusion every 3 weeks:",
      "0.03-3.4 mg/kg in DREAMM-1 (1-hour infusion), 2.5 or 3.4 mg/kg in",
      "DREAMM-2 (30-minute infusion), 2.5 mg/kg in DREAMM-3 and the",
      "DREAMM-5 monotherapy arm (infusion over at least 30 minutes).",
      sep = " "
    ),
    regions = "Multinational",
    renal_function = "Normal 31.2%, mild impairment 42.7%, moderate 23.7%, severe 2.2% (Table 1).",
    biomarker_summary = paste(
      "Baseline medians (ranges), Table 1 overall: sBCMA 72.9 (2.08-2030)",
      "ng/mL; M-protein 18.0 (0.70-102) g/L; albumin 38.0 (19.0-52.3) g/L;",
      "IgG 12.0 (0.30-119) g/L; beta-2-microglobulin 320 (102-1970) nmol/L.",
      sep = " "
    ),
    notes = paste(
      "Table 1: DREAMM-1 (NCT02064387) n = 45, DREAMM-2 (NCT03525678)",
      "n = 216, DREAMM-3 (NCT04162210) n = 215, DREAMM-5 (NCT04126200)",
      "monotherapy arm n = 34. 15,860 observations: 4586 total ADC, 4624",
      "free sBCMA, 3422 complexed sBCMA and 3228 M-protein. Myeloma",
      "immunoglobulin type IgG 75.5%, IgA 21.6%.",
      sep = " "
    ),
    software = "NONMEM 7.4.3 (SAEM followed by importance sampling)."
  )

  ini({
    # ------------------------------------------------------------------
    # Structural parameters -- Clements 2026 Table 2 (final model).
    # Reference subject: male, 74.65 kg, albumin 38 g/L, sBCMA 72.85 ng/mL,
    # no extramedullary disease.
    # ------------------------------------------------------------------
    lcl <- log(0.7); label("Clearance of total belantamab mafodotin, CLADC (L/day)") # Table 2: CLADC = 0.7 L/day (RSE 5.12%)
    lvc <- log(4.81); label("Central volume of distribution of total belantamab mafodotin, VC (L)") # Table 2: VC = 4.81 L (RSE 1.83%)
    lq <- log(0.771); label("Intercompartmental clearance of total belantamab mafodotin, QADC (L/day)") # Table 2: QADC = 0.771 L/day (RSE 6.98%)
    lvp <- log(12.7); label("Peripheral volume of distribution of total belantamab mafodotin, VP (L)") # Table 2: VP = 12.7 L (RSE 9.06%)
    lkdeg <- log(0.671); label("Free sBCMA degradation rate constant, KDEG (1/day)") # Table 2: KDEG = 0.671 1/day (RSE 9.69%)
    lkd <- log(1.11); label("Belantamab mafodotin-sBCMA equilibrium dissociation constant, KD (nM)") # Table 2: KD = 1.11 'nM/L' (RSE 3.85%)
    lkkill <- log(0.00019); label("Second-order M-protein kill rate constant per unit free drug, KDEATH (1/(nM*day))") # Table 2: KDEATH = 0.00019 1/(nM*day) (RSE 10.3%)
    lp <- log(0.00319); label("M-protein logistic growth rate constant, KGR (1/day)") # Table 2: KGR = 0.00319 1/day (RSE 17.2%)
    lkpout <- log(4.45); label("Rate constant for sBCMA diffusion out of the reservoir, KDIFF (1/day)") # Table 2: KDIFF = 4.45 1/day (RSE 9.55%)

    # ------------------------------------------------------------------
    # Covariate effects -- Clements 2026 Table 2. Continuous covariates
    # enter as (X / median)^theta and categorical ones as (1 + theta)
    # (Supporting Information, 'Covariate mathematical structures';
    # control stream $PK).
    # ------------------------------------------------------------------
    e_alb_cl <- -0.879; label("Power exponent of baseline albumin on CL, ALB/38 (unitless)") # Table 2: Baseline albumin effect on CLADC = -0.879 (RSE 17.8%)
    e_sbcma_cl <- 0.258; label("Power exponent of baseline sBCMA on CL, SBCMA/72.85 (unitless)") # Table 2: Baseline sBCMA effect on CLADC = 0.258 (RSE 13.5%)
    e_sexf_cl <- -0.126; label("Proportional shift in CL for females (fraction)") # Table 2: Effect of female sex on CLADC = -0.126 (RSE 32.1%)
    e_wt_cl <- 0.488; label("Power exponent of baseline body weight on CL, WT/74.65 (unitless)") # Table 2: Body weight effect on CLADC = 0.488 (RSE 19.7%)
    e_sbcma_kkill <- 0.318; label("Power exponent of baseline sBCMA on KDEATH, SBCMA/72.85 (unitless)") # Table 2: Baseline sBCMA effect on KDEATH = 0.318 (RSE 19.6%)
    e_sexf_vc <- -0.136; label("Proportional shift in VC for females (fraction)") # Table 2: Effect of female sex on belantamab mafodotin VC = -0.136 (RSE 19.4%)
    e_wt_vc <- 0.506; label("Power exponent of baseline body weight on VC, WT/74.65 (unitless)") # Table 2: Body weight effect on belantamab mafodotin VC = 0.506 (RSE 13.1%)
    e_dis_emd_kkill <- -0.684; label("Proportional shift in KDEATH with extramedullary disease (fraction)") # Table 2: Extramedullary disease effect on KDEATH = -0.684 (RSE 12.7%)

    # ------------------------------------------------------------------
    # Between-subject variability. Table 2 reports CV% = 100 *
    # sqrt(exp(omega^2) - 1) (each row's SE / RSE reproduces this
    # variance), so omega^2 = log(1 + CV^2).
    # ------------------------------------------------------------------
    etalcl ~ 0.15188 # Table 2: BSV CV% CLADC = 40.5 -> log(1 + 0.405^2)
    etalvc ~ 0.06300 # Table 2: BSV CV% VC = 25.5 -> log(1 + 0.255^2)
    etalvp ~ 1.03762 # Table 2: BSV CV% VP = 135 -> log(1 + 1.35^2)
    etalkkill ~ 1.02805 # Table 2: BSV CV% KDEATH = 134 -> log(1 + 1.34^2)
    etalp ~ 2.38505 # Table 2: BSV CV% KGR = 314 -> log(1 + 3.14^2)
    # The control stream also carries etas on QA, KDEG, KD and KDIFF with
    # variances fixed to 0.0225 (15% CV); Table 2 lists only estimated
    # variances, so these come from the Supporting Information $OMEGA.
    etalq ~ fixed(0.0225) # Supporting Information control stream $OMEGA: 0.0225 FIX ; QA
    etalkdeg ~ fixed(0.0225) # Supporting Information control stream $OMEGA: 0.0225 FIX ; KDEG
    etalkd ~ fixed(0.0225) # Supporting Information control stream $OMEGA: 0.0225 FIX ; KD
    etalkpout ~ fixed(0.0225) # Supporting Information control stream $OMEGA: 0.0225 FIX ; KDIFF

    # ------------------------------------------------------------------
    # Residual error: additive on the natural-log scale for every analyte
    # (Table 2 reports the variance; lnorm() takes the SD).
    # ------------------------------------------------------------------
    expSd_sbcma_free <- 0.65115; label("Free sBCMA residual SD on the log scale (unitless)") # Table 2: residual variance free sBCMA = 0.424 -> sqrt
    expSd_sbcma_complex <- 0.78038; label("Complexed sBCMA residual SD on the log scale (unitless)") # Table 2: residual variance complex sBCMA = 0.609 -> sqrt
    expSd <- 0.36194; label("Total belantamab mafodotin residual SD on the log scale (unitless)") # Table 2: residual variance total belantamab mafodotin = 0.131 -> sqrt
    expSd_mprotein_obs <- 0.32863; label("Serum M-protein residual SD on the log scale (unitless)") # Table 2: residual variance serum M-protein = 0.108 -> sqrt
  })

  model({
    # Molecular weights (kDa) from the control stream: MWADC = 152.1,
    # MWSBCMA = 6. mg/L * 1000 / MW(kDa) = nM; ng/mL / MW(kDa) = nM.
    mw_adc <- 152.1
    mw_sbcma <- 6

    # 1. Individual parameters (control stream $PK)
    cl <- exp(lcl + etalcl) * (ALB / 38)^e_alb_cl * (SBCMA / 72.85)^e_sbcma_cl *
      (1 + e_sexf_cl * SEXF) * (WT / 74.65)^e_wt_cl
    vc <- exp(lvc + etalvc) * (1 + e_sexf_vc * SEXF) * (WT / 74.65)^e_wt_vc
    q <- exp(lq + etalq)
    vp <- exp(lvp + etalvp)
    kdeg <- exp(lkdeg + etalkdeg)
    kd <- exp(lkd + etalkd)
    kkill <- exp(lkkill + etalkkill) * (SBCMA / 72.85)^e_sbcma_kkill *
      (1 + e_dis_emd_kkill * DIS_EMD)
    p <- exp(lp + etalp)
    kpout <- exp(lkpout + etalkpout)

    # 2. Micro-constants. KINT is forced equal to KEL (Results; Figure 2).
    kel <- cl / vc
    kint <- kel
    k12 <- q / vc
    k21 <- q / vp

    # Baseline sBCMA in nM; synthesis rate KSYN = KDEG * SBCMABL / MWSBCMA.
    sbcma0 <- SBCMA / mw_sbcma
    ksyn <- kdeg * sbcma0

    # 3. Quasi-equilibrium binding (control stream $DES). ctot is total
    # belantamab mafodotin and cfree the free drug, both nM; total_target
    # is total sBCMA in nM.
    ctot <- central / vc * 1000 / mw_adc
    bt <- ctot - total_target - kd
    cfree <- 0.5 * (bt + sqrt(bt^2 + 4 * kd * ctot))
    cfree <- max(cfree, 0)
    mp <- max(mprotein, 0)

    # 4. ODEs. The nmol/day binding fluxes are converted back to mg/day for
    # the mg-valued ADC state. Peripheral distribution acts on total drug.
    d/dt(central) <- -(kel * cfree * vc + kint * total_target * cfree / (kd + cfree) * vc) * mw_adc / 1000 -
      k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(total_target) <- kpout * precursor1 - (kint - kdeg) * total_target * cfree / (kd + cfree) -
      kdeg * total_target
    # Logistic M-protein growth with the 150 g/L carrying capacity hard-coded
    # in the control stream, and a kill term linear in free drug.
    d/dt(mprotein) <- p * (1 - mp / 150) * mp - kkill * cfree * mp
    d/dt(precursor1) <- ksyn * (mp / MCPROT) - kpout * precursor1

    # Initial conditions (control stream A_0). The reservoir starts at the
    # baseline sBCMA concentration (Methods, 'Model Finalization'), not at
    # its drug-free steady state KSYN / KDIFF.
    total_target(0) <- sbcma0
    mprotein(0) <- MCPROT
    precursor1(0) <- sbcma0

    # 5. Observation model (control stream $ERROR). Free and complexed
    # drug and sBCMA are recomputed from the states; free drug is not
    # floored here, matching the control stream.
    cf_obs <- 0.5 * (bt + sqrt(bt^2 + 4 * kd * ctot))
    sbcma_cplx_true <- total_target * cf_obs / (kd + cf_obs)
    sbcma_free_true <- total_target - sbcma_cplx_true

    # Assay-dilution correction: on dilution by a factor DIL the complex
    # partly re-equilibrates, a fraction alpha dissociating. Complexed
    # sBCMA is measured at DIL = 10 and free sBCMA at DIL = 20.
    if (sbcma_cplx_true == 0) {
      alpha10 <- 0
      alpha20 <- 0
    } else {
      bq10 <- sbcma_free_true + cf_obs + 10 * kd
      alpha10 <- (-bq10 + sqrt(bq10^2 + 4 * sbcma_cplx_true * (10 - 1) * kd)) / (2 * sbcma_cplx_true)
      bq20 <- sbcma_free_true + cf_obs + 20 * kd
      alpha20 <- (-bq20 + sqrt(bq20^2 + 4 * sbcma_cplx_true * (20 - 1) * kd)) / (2 * sbcma_cplx_true)
    }
    sbcma_complex_dil <- sbcma_cplx_true * (1 - alpha10)
    sbcma_free_dil <- sbcma_free_true + alpha20 * sbcma_cplx_true

    # Lower floors applied before the log transform (control stream):
    # free sBCMA 0.5 nM, complexed sBCMA 0.1 nM, total ADC 1 nM,
    # M-protein 1 g/L.
    sbcma_free <- max(sbcma_free_dil, 0.5)
    sbcma_complex <- max(sbcma_complex_dil, 0.1)
    Cc <- max(ctot, 1)
    mprotein_obs <- max(mp, 1)

    sbcma_free ~ lnorm(expSd_sbcma_free)
    sbcma_complex ~ lnorm(expSd_sbcma_complex)
    Cc ~ lnorm(expSd)
    mprotein_obs ~ lnorm(expSd_mprotein_obs)
  })
}
