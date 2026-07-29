Williams_2016_rituximab_das28cfb <- function() {
  description <- paste(
    "Longitudinal PD model of Disease Activity Score in 28 joints change from",
    "baseline (DAS28cfb) for rituximab (both the proposed biosimilar PF-05280586",
    "and the two rituximab reference products sourced from the EU and US) in",
    "adults with rheumatoid arthritis on background methotrexate with prior",
    "inadequate response to one or more TNF-antagonist therapies (Williams 2016).",
    "The model is an inhibitory-Emax-on-log-scale exposure-response for the",
    "concentration-driven drug effect superimposed on an exponential-onset",
    "placebo/background time course, with the outer transformation",
    "DAS28cfb = 1 - exp(fnon-C(t) + fC(C)) so that improvement is a negative",
    "change from baseline. The typical maximum placebo/background effect PMAX",
    "and typical maximum drug effect Emax are log-scale additive parameters",
    "with correlated between-subject etas; a separate additive subject-level",
    "eta on the outer scale (etaCS) captures the within-subject correlation",
    "induced by computing change from baseline. Covariate model includes the",
    "individual components of baseline disease activity (TJ28, SJ28, log CRP+1,",
    "patient global assessment) and treatment-arm indicators on each of the",
    "three main structural parameters (PMAX, kp, Emax) with rituximab-EU as",
    "the reference arm. This is a PD-only extraction: the rituximab plasma",
    "concentration Cc is supplied as a time-varying covariate CP_RITUXIMAB_UGML",
    "(the Williams 2016 popPK model - the source of the individual predicted",
    "concentrations Cij feeding the DAS28cfb model - was a two-compartment",
    "structural model whose parameter values are not reported in the paper or",
    "its supplements; users wanting to drive this PD model from a simulated PK",
    "source must supply their own concentration trajectory).",
    sep = " "
  )

  reference <- paste(
    "Williams JH, Hutmacher MM, Zierhut ML, Becker JC, Gumbiner B,",
    "Spencer-Green G, Melia LA, Liao KH, Suster M, Yin D, Li R, Meng X.",
    "Comparative assessment of clinical response in patients with rheumatoid",
    "arthritis between PF-05280586, a proposed rituximab biosimilar, and",
    "rituximab. Br J Clin Pharmacol. 2016;82(6):1568-1579.",
    "doi:10.1111/bcp.13094. PMID: 27530379.",
    "Companion ACR responder rate model from the same paper:",
    "modellib('Williams_2016_rituximab_acr').",
    "The main paper's popPK model was described as a two-compartment structural",
    "model with baseline body surface area and sex on clearance and central",
    "volume, similar to the earlier rituximab popPK of",
    "Ng CM, Bruno R, Combs D, Davies B.",
    "Population pharmacokinetics of rituximab (anti-CD20 monoclonal antibody)",
    "in rheumatoid arthritis patients during a phase II clinical trial.",
    "J Clin Pharmacol 2005;45(7):792-801. doi:10.1177/0091270005277075.",
    "Williams 2016 does not tabulate the popPK parameter estimates - the",
    "in-model concentration driver Cij is a paper-internal individual",
    "posthoc from that unreported popPK fit.",
    sep = " "
  )

  vignette <- "Williams_2016_rituximab_biosimilar"

  paper_specific_etas <- c("etaCS")
  paper_specific_residual_sds <- c("addSd_das28")

  units <- list(
    time          = "week",
    dosing        = "(none; PD-only model fed by an external rituximab plasma-concentration covariate)",
    concentration = "(observation DAS28cfb is the DAS28-CRP change from baseline in DAS28 units; driving covariate CP_RITUXIMAB_UGML is in ug/mL)"
  )

  covariateData <- list(
    CP_RITUXIMAB_UGML = list(
      description        = "Instantaneous rituximab (or biosimilar PF-05280586) plasma concentration at each PD observation, supplied as a time-varying covariate from observed serum samples or an upstream popPK source.",
      units              = "ug/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying per event row. Drives the concentration-effect Emax term",
        "fC(C) = Emax_i * C / (C + exp(lec50)) inside model(). Williams 2016",
        "used individual predicted concentrations Cij from a two-compartment",
        "popPK model that included baseline body-surface-area and sex effects",
        "on CL and Vc (Methods 'Population PK/PD models' section); the popPK",
        "parameter values are not tabulated in the paper or supplements and",
        "must be supplied externally. Typical steady-state Cmax after the",
        "two-dose loading course (1000 mg IV on day 1 and day 15) reaches the",
        "hundreds of ug/mL and washes out to <1 ug/mL by week 24 - see",
        "Williams 2016 Supplemental Figure S1 (VPC of the rituximab popPK",
        "model). Set to 0 outside the drug-exposure window; the concentration",
        "Emax term then collapses to 0.",
        sep = " "
      ),
      source_name        = "CONC"
    ),
    TEND_28JOINT = list(
      description        = "Baseline tender joint count on the 28-joint DAS28 scale (integer 0-28).",
      units              = "count (0-28)",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject. Enters the covariate model as an additive",
        "log-scale effect (TEND_28JOINT - 16) on each of PMAX, kp, and Emax",
        "(Williams 2016 Supplemental Methods; reference value 16 = paper's",
        "declared median of the DAS28cfb dataset). Component of the DAS28",
        "composite score.",
        sep = " "
      ),
      source_name        = "TJ28"
    ),
    SWOL_28JOINT = list(
      description        = "Baseline swollen joint count on the 28-joint DAS28 scale (integer 0-28).",
      units              = "count (0-28)",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject. Enters the covariate model as an additive",
        "log-scale effect (SWOL_28JOINT - 12) on each of PMAX, kp, and Emax",
        "(Williams 2016 Supplemental Methods; reference value 12 = paper's",
        "declared median of the DAS28cfb dataset).",
        sep = " "
      ),
      source_name        = "SJ28"
    ),
    CRP = list(
      description        = "Baseline C-reactive protein (BCRP).",
      units              = "mg/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject (baseline value). Enters the covariate model",
        "on the log(CRP + 1) scale to accommodate the highly skewed",
        "distribution and BCRP = 0 observations (Williams 2016 Supplemental",
        "Methods): lBCRP_i = log(CRP + 1). The reference value on that",
        "log-transformed scale is 2.2 log-units (= log(9.03 - 1 + 1)",
        "approximately, matching the paper's declared median lBCRP). Applied",
        "as (log(CRP + 1) - 2.2) additively on PMAX, kp, and Emax.",
        sep = " "
      ),
      source_name        = "BCRP"
    ),
    PGA_PT = list(
      description        = "Baseline patient's global assessment of arthritis (100-mm visual analogue scale).",
      units              = "mm (0-100 VAS)",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject. Additive log-scale effect (PGA_PT - 70) on",
        "each of PMAX, kp, and Emax; reference 70 mm is the paper's declared",
        "median of the DAS28cfb dataset (Williams 2016 Supplemental Methods).",
        "Distinct from the physician's global assessment (BLPHYVAS) and from",
        "patient-reported pain (PAIN); PGA_PT is the patient's own overall",
        "rating of arthritis disease activity.",
        sep = " "
      ),
      source_name        = "PGA"
    ),
    TRT = list(
      description        = "Treatment-arm integer indicator: 0 = rituximab-EU (reference), 1 = PF-05280586 (proposed biosimilar), 2 = rituximab-US.",
      units              = "(categorical / integer-coded)",
      type               = "categorical",
      reference_category = "0 (rituximab-EU)",
      notes              = paste(
        "Time-fixed per subject. Encodes the trial's three treatment arms.",
        "The paper's covariate model estimates one indicator effect per arm",
        "on each of PMAX, kp, and Emax with rituximab-EU as the reference",
        "level (Williams 2016 Supplemental Methods). Two derived indicators",
        "are computed inside model() as TRT_PF = 1*(TRT == 1) and",
        "TRT_US = 1*(TRT == 2), each entering the covariate function",
        "additively on the log scale.",
        sep = " "
      ),
      source_name        = "TRT"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 214L,
    n_studies        = 1L,
    n_observations   = 1382L,
    age_range        = "adults (>=18 years); mean (SD) 54.8 (11.7) years (PF-05280586), 55.7 (10.2) years (rituximab-EU), 53.8 (11.8) years (rituximab-US)",
    weight_range     = "mean (SD) 86.2 (22.0) kg (PF-05280586), 82.6 (19.8) kg (rituximab-EU), 80.4 (21.6) kg (rituximab-US)",
    sex_female_pct   = 77.6,
    disease_state    = "Active rheumatoid arthritis on background methotrexate with inadequate response to one or more TNF-antagonist therapies. Baseline DAS28-CRP mean (SD): 5.64 (0.85) PF-05280586, 5.80 (0.96) rituximab-EU, 6.2 (0.89) rituximab-US.",
    dose_range       = "1000 mg IV on days 1 and 15 (standard rituximab RA induction course). All subjects received 100 mg IV methylprednisolone premedication.",
    regions          = "Multi-regional biosimilar development trial (ClinicalTrials.gov NCT01526057).",
    baseline_disease = list(
      DAS28_CRP_by_arm = "PF-05280586 5.64 (0.85); rituximab-EU 5.80 (0.96); rituximab-US 6.2 (0.89)",
      SWOL_28JOINT_by_arm = "PF-05280586 11.4 (5.0); rituximab-EU 13.0 (6.6); rituximab-US 14.0 (6.0)",
      TEND_28JOINT_by_arm = "PF-05280586 14.3 (6.5); rituximab-EU 15.1 (6.8); rituximab-US 18.0 (6.5)",
      CRP_by_arm_mg_L    = "PF-05280586 12.4 (14.9); rituximab-EU 14.7 (17.6); rituximab-US 18.2 (25.1)",
      PGA_by_arm         = "PF-05280586 67.4 (16.8); rituximab-EU 67.7 (20.9); rituximab-US 74.8 (16.0)"
    ),
    disease_duration = "Mean (SD) since first RA diagnosis: 12.7 (8.4), 11.8 (8.3), 10.6 (8.1) years for PF-05280586, rituximab-EU, rituximab-US arms respectively.",
    notes            = paste(
      "Baseline demographics from Williams 2016 Table 1. The DAS28cfb dataset",
      "included 214 baseline and 1382 postbaseline observations",
      "(Williams 2016 Results 'Population PK/PD models'). Rituximab-US arm",
      "had the highest baseline DAS28 with numerically more severe disease",
      "and shorter duration of RA - baseline imbalance across arms was the",
      "motivation for including the individual DAS28 components as covariates",
      "rather than pooling. Ethnicity, race, and country/region were not",
      "reported quantitatively in Table 1 for the popPK/PD analysis population.",
      sep = " "
    )
  )

  ini({
    # --------------------------------------------------------------------------
    # Structural parameters (Williams 2016 Supplemental Table S1, With
    # Baseline Covariates column). All THETA values are on the log scale by
    # design so that covariate effects and etas enter additively; the equation
    # DAS28cfb = 1 - exp(fnon-C + fC) exponentiates back to the DAS28 scale
    # so PMAX and Emax function as multiplicative modifiers of the outer
    # 1 - exp() transformation. Transformed estimates given in parentheses
    # for reader convenience.
    # --------------------------------------------------------------------------
    lpmax     <-  1.08;    label("Typical PMAX for placebo/background exponential-onset effect (log DAS28-scale; typical PMAX = exp(1.08) = 2.95)") # Williams 2016 Supp Table S1: PMAX estimate 1.08 (SE 0.0603)
    lemax     <- -0.0579;  label("Typical Emax for concentration-driven drug effect (log DAS28-scale fractional; typical Emax = exp(-0.0579) = 0.944)") # Williams 2016 Supp Table S1: Emax estimate -0.0579 (SE 0.107)
    lthalfrec <-  0.466;   label("Log typical onset half-life for the placebo/background exponential rise, giving kp = log(2)/exp(lthalfrec) (weeks)") # Williams 2016 Supp Table S1: Onset (kp) estimate 0.466 (SE 0.238); Onset labelled as half-life in weeks
    lec50     <-  2.20;    label("Log typical rituximab concentration at half-maximal Emax (ug/mL)") # Williams 2016 Supp Table S1: EC50 estimate 2.20 (SE 0.798); transformed = exp(2.20) = 9.05 ug/mL

    # --------------------------------------------------------------------------
    # Covariate effects on PMAX (Williams 2016 Supp Table S1; all additive on
    # log-DAS28 scale; interpretation: exp(e) is the multiplicative fold-change
    # of the transformed PMAX per unit covariate deviation from reference).
    # --------------------------------------------------------------------------
    e_tj28_pmax    <-  0.00347;   label("Log-scale additive effect of TEND_28JOINT - 16 on PMAX (per count)")  # Williams 2016 Supp Table S1: TJ28 on PMAX 0.00347 (SE 0.00756)
    e_sj28_pmax    <-  0.00619;   label("Log-scale additive effect of SWOL_28JOINT - 12 on PMAX (per count)")   # Williams 2016 Supp Table S1: SJ28 on PMAX 0.00619 (SE 0.00872)
    e_lbcrp_pmax   <- -0.00944;   label("Log-scale additive effect of (log(CRP+1) - 2.2) on PMAX (per log unit)") # Williams 2016 Supp Table S1: lBCRP on PMAX -0.00944 (SE 0.0458)
    e_pga_pmax     <-  0.00424;   label("Log-scale additive effect of PGA_PT - 70 on PMAX (per mm VAS)")         # Williams 2016 Supp Table S1: PGA on PMAX 0.00424 (SE 0.00165)
    e_trt_pf_pmax  <- -0.0951;    label("Log-scale additive effect on PMAX for TRT = 1 (PF-05280586)")           # Williams 2016 Supp Table S1: PF-05280586 on PMAX -0.0951 (SE 0.0835)
    e_trt_us_pmax  <-  0.0938;    label("Log-scale additive effect on PMAX for TRT = 2 (rituximab-US)")          # Williams 2016 Supp Table S1: Ritux-US on PMAX 0.0938 (SE 0.0772)

    # --------------------------------------------------------------------------
    # Covariate effects on the onset half-life (Williams 2016 Supp Table S1;
    # positive value means longer half-life = slower onset).
    # --------------------------------------------------------------------------
    e_tj28_thalfrec    <-  0.00365;  label("Log-scale additive effect of TEND_28JOINT - 16 on lthalfrec")  # Williams 2016 Supp Table S1: TJ28 on kp 0.00365 (SE 0.0354)
    e_sj28_thalfrec    <-  0.00900;  label("Log-scale additive effect of SWOL_28JOINT - 12 on lthalfrec")   # Williams 2016 Supp Table S1: SJ28 on kp 0.00900 (SE 0.0283)
    e_lbcrp_thalfrec   <-  0.266;    label("Log-scale additive effect of (log(CRP+1) - 2.2) on lthalfrec") # Williams 2016 Supp Table S1: lBCRP on kp 0.266 (SE 0.180)
    e_pga_thalfrec     <- -0.00452;  label("Log-scale additive effect of PGA_PT - 70 on lthalfrec")         # Williams 2016 Supp Table S1: PGA on kp -0.00452 (SE 0.00597)
    e_trt_pf_thalfrec  <-  0.123;    label("Log-scale additive effect on lthalfrec for TRT = 1 (PF-05280586)") # Williams 2016 Supp Table S1: PF-05280586 on kp 0.123 (SE 0.300)
    e_trt_us_thalfrec  <- -0.144;    label("Log-scale additive effect on lthalfrec for TRT = 2 (rituximab-US)") # Williams 2016 Supp Table S1: Ritux-US on kp -0.144 (SE 0.300)

    # --------------------------------------------------------------------------
    # Covariate effects on Emax (Williams 2016 Supp Table S1).
    # --------------------------------------------------------------------------
    e_tj28_emax    <- -0.00560;   label("Log-scale additive effect of TEND_28JOINT - 16 on Emax (per count)")   # Williams 2016 Supp Table S1: TJ28 on Emax -0.00560 (SE 0.0115)
    e_sj28_emax    <-  0.00507;   label("Log-scale additive effect of SWOL_28JOINT - 12 on Emax (per count)")    # Williams 2016 Supp Table S1: SJ28 on Emax 0.00507 (SE 0.0126)
    e_lbcrp_emax   <-  0.0839;    label("Log-scale additive effect of (log(CRP+1) - 2.2) on Emax (per log unit)") # Williams 2016 Supp Table S1: lBCRP on Emax 0.0839 (SE 0.0698)
    e_pga_emax     <- -0.00186;   label("Log-scale additive effect of PGA_PT - 70 on Emax (per mm VAS)")         # Williams 2016 Supp Table S1: PGA on Emax -0.00186 (SE 0.00306)
    e_trt_pf_emax  <-  0.133;     label("Log-scale additive effect on Emax for TRT = 1 (PF-05280586)")           # Williams 2016 Supp Table S1: PF-05280586 on Emax 0.133 (SE 0.131)
    e_trt_us_emax  <- -0.0350;    label("Log-scale additive effect on Emax for TRT = 2 (rituximab-US)")          # Williams 2016 Supp Table S1: Ritux-US on Emax -0.0350 (SE 0.145)

    # --------------------------------------------------------------------------
    # Inter-individual variability on PMAX and Emax (correlated block).
    # Williams 2016 Supp Table S1 reports:
    #   Var(eta1 for PMAX)  transformed 32.5% CV -> omega^2 = log(1 + 0.325^2) = 0.1004
    #   Var(eta2 for Emax)  transformed 31.4% CV -> omega^2 = log(1 + 0.314^2) = 0.0941
    #   rho(eta1, eta2)     transformed -0.301   -> cov = -0.301 * sqrt(0.1004 * 0.0941) = -0.0292
    # kp has no eta (paper equation: kp = log(2)/exp(lthalfrec + g_i(theta))).
    # --------------------------------------------------------------------------
    etalpmax + etalemax ~ c(0.1004, -0.0292, 0.0941)

    # --------------------------------------------------------------------------
    # Subject-level additive random effect on the outer 1 - exp() equation
    # (Williams 2016 Supp Methods: 'the correlation in observations at the
    # residual level induced by computing change from baseline'). Reported as
    # SD 0.475 DAS28 units on the transformed scale, so variance = 0.475^2 =
    # 0.2256 DAS28^2. Paper-specific eta (not tied to a log-transformed
    # structural parameter).
    # --------------------------------------------------------------------------
    etaCS ~ 0.2256

    # --------------------------------------------------------------------------
    # Additive residual error on DAS28cfb (Williams 2016 Supp Table S1
    # transformed sigma = 0.693 DAS28 units).
    # --------------------------------------------------------------------------
    addSd_das28 <- 0.693;   label("Additive residual SD on DAS28cfb (DAS28 units)") # Williams 2016 Supp Table S1: sigma estimate -0.366 (log SD); transformed 0.693 DAS28
  })

  model({
    # Treatment-arm indicators (per Williams 2016 Supp Methods
    # 'g_i^(A)(theta) = ... + theta_PF I(TRT_i = PF-05280586) + theta_US I(TRT_i = Ritux-US)').
    TRT_PF <- (TRT == 1)
    TRT_US <- (TRT == 2)

    # Log-transformed baseline CRP (Williams 2016 Supp Methods:
    # lBCRP = log(BCRP + 1); reference median 2.2 log-units).
    lBCRP <- log(CRP + 1)

    # Individual log-scale structural parameters with covariate effects
    # centred on the paper's declared reference values (TJ28 = 16, SJ28 = 12,
    # lBCRP = 2.2, PGA = 70, TRT = 0 = rituximab-EU).
    pmax_i <- lpmax + etalpmax +
              e_tj28_pmax    * (TEND_28JOINT - 16) +
              e_sj28_pmax    * (SWOL_28JOINT - 12) +
              e_lbcrp_pmax   * (lBCRP - 2.2) +
              e_pga_pmax     * (PGA_PT - 70) +
              e_trt_pf_pmax  * TRT_PF +
              e_trt_us_pmax  * TRT_US

    lthalfrec_i <- lthalfrec +
                    e_tj28_thalfrec    * (TEND_28JOINT - 16) +
                    e_sj28_thalfrec    * (SWOL_28JOINT - 12) +
                    e_lbcrp_thalfrec   * (lBCRP - 2.2) +
                    e_pga_thalfrec     * (PGA_PT - 70) +
                    e_trt_pf_thalfrec  * TRT_PF +
                    e_trt_us_thalfrec  * TRT_US
    kp_i <- log(2) / exp(lthalfrec_i)

    emax_i <- lemax + etalemax +
              e_tj28_emax    * (TEND_28JOINT - 16) +
              e_sj28_emax    * (SWOL_28JOINT - 12) +
              e_lbcrp_emax   * (lBCRP - 2.2) +
              e_pga_emax     * (PGA_PT - 70) +
              e_trt_pf_emax  * TRT_PF +
              e_trt_us_emax  * TRT_US

    ec50_i <- exp(lec50)

    # Model equations (Williams 2016 Supp Methods 'DAS28cfb model'):
    #   fnon-C(t) = PMAX_i * (1 - exp(-kp_i * t))
    #   fC(C)     = Emax_i * C / (C + exp(lec50))
    #   DAS28cfb  = 1 - exp(fnon-C + fC) + etaCS + eps
    # PMAX and Emax enter the exponent additively on the log DAS28 scale;
    # exponentiating gives the multiplicative transformed values shown in
    # Supp Table S1 (typical PMAX = 2.95 DAS28, typical Emax = 0.944 fractional).
    fnonC <- pmax_i * (1 - exp(-kp_i * time))
    conc  <- CP_RITUXIMAB_UGML
    fC    <- emax_i * conc / (conc + ec50_i)
    das28cfb <- 1 - exp(fnonC + fC) + etaCS

    das28cfb ~ add(addSd_das28)
  })
}
