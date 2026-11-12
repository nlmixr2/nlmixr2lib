Hodson_2023_radiotherapy_antipdl1_olaparib_azd0156_mouse <- function() {
  description <- paste(
    "Preclinical (mouse, female C57BL/6 with subcutaneous MC38 tumors). QSP.",
    "Tumor-immune model of radiotherapy (RT) combined with anti-PD-L1 and/or",
    "a DNA-damage-response inhibitor (olaparib or AZD0156): a logistic",
    "proliferating rim and a quiescent core make up the tumor diameter, RT",
    "recruits antigen-presenting cells that activate cytolytic CD8 T cells,",
    "and each combination drug lowers the T-cell exhaustion rate."
  )
  reference <- paste(
    "Hodson D, Mistry H, Guzzetti S, Davies M, Staniszewska A, Farrington P,",
    "Cadogan E, Yates J, Aarons L, Ogungbenro K. Mixed effects modeling of radiotherapy in combination",
    "with immune checkpoint blockade or inhibitors of the DNA damage",
    "response pathway. CPT Pharmacometrics Syst Pharmacol.",
    "2023;12(11):1640-1652. doi:10.1002/psp4.13026. PMID 37722071;",
    "PMCID PMC10681475."
  )
  vignette <- "Hodson_2023_radiotherapy_antipdl1_olaparib_azd0156_mouse"
  # rim / core hold tumor depths (mm); apc_cells / cd8_cells are the paper's
  # A(t) and T(t) immune pools. None generalises beyond this model.
  paper_specific_compartments <- c("rim", "core", "apc_cells", "cd8_cells")
  units <- list(
    time = "day",
    dosing = "none (treatments enter as TRT_* indicators; RT as a Heaviside input in Gy)",
    concentration = "mm (tumor diameter; not a drug concentration)"
  )

  compartmentData <- list(
    rim = list(analyte = "proliferating tumor rim depth", units = "mm", specimen = "not applicable", verified = TRUE),
    core = list(analyte = "quiescent tumor core depth", units = "mm", specimen = "not applicable", verified = TRUE),
    apc_cells = list(
      analyte = "antigen-presenting cells",
      units = "cells/mm",
      specimen = "not applicable",
      verified = TRUE
    ),
    cd8_cells = list(
      analyte = "active cytolytic CD8 T cells",
      units = "cells/mm",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  covariateData <- list(
    TRT_RADIOTHERAPY = list(
      description = "Radiotherapy group indicator (1 = mouse received 2 Gy x 5 daily fractions from day 0)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no radiotherapy; control cohort)",
      notes = paste(
        "Time-fixed per mouse. Switches on the Heaviside RT input",
        "R(t) = rt_dose * (t < rt_end) of Equation 2. In the NONMEM code",
        "(Code S1) this is the flag RT1 = 1 for TRT > 1; every fitted cohort",
        "other than the control received RT."
      ),
      source_name = "RT1 (derived from TRT)"
    ),
    TRT_ANTIPDL1 = list(
      description = "Anti-mouse-PD-L1 (ICI) group indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no anti-PD-L1)",
      notes = paste(
        "Time-fixed per mouse; 1 for the RT/ICI, RT/PARPi/ICI and RT/ATMi/ICI",
        "cohorts (10 mg/kg IP twice weekly for 3 weeks). Divides the T-cell",
        "exhaustion rate by (1 + gamma_ICI + ...). Code S1 flag IO1."
      ),
      source_name = "IO1 (derived from TRT)"
    ),
    TRT_OLAPARIB = list(
      description = "Olaparib (PARP inhibitor) group indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no olaparib)",
      notes = paste(
        "Time-fixed per mouse; 1 for the RT/PARPi and RT/PARPi/ICI cohorts",
        "(100 mg/kg PO daily for 21 days). Code S1 flag PARP."
      ),
      source_name = "PARP (derived from TRT)"
    ),
    TRT_AZD0156 = list(
      description = "AZD0156 (ATM inhibitor) group indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no AZD0156)",
      notes = paste(
        "Time-fixed per mouse; 1 for the RT/ATMi and RT/ATMi/ICI cohorts",
        "(2.5 mg/kg PO daily for 21 days). Code S1 flag ATM1. Mutually",
        "exclusive with TRT_OLAPARIB in the source studies."
      ),
      source_name = "ATM1 (derived from TRT)"
    )
  )

  population <- list(
    species = "mouse (female C57BL/6; subcutaneous MC38 syngeneic colorectal tumor)",
    n_subjects = 212L,
    n_studies = 4L,
    sex_female_pct = 100,
    disease_state = "Subcutaneous MC38 murine colorectal carcinoma (10^7 cells injected 3 days before treatment start)",
    dose_range = paste(
      "RT 2 Gy/day on days 0-4 (10 Gy total); anti-mouse-PD-L1 10 mg/kg IP",
      "twice weekly for 3 weeks; olaparib 100 mg/kg PO daily for 21 days;",
      "AZD0156 2.5 mg/kg PO daily for 21 days."
    ),
    regions = "Preclinical (UK; AstraZeneca DDRIO studies)",
    notes = paste(
      "Four in vivo experiments (DDRIO1815, DDRIO1821: RT +/- olaparib +/-",
      "anti-PD-L1; DDRIO1824, DDRIO1835: RT +/- AZD0156 +/- anti-PD-L1;",
      "Table 1). The model was fit only to the control, RT, RT/DDRi, RT/ICI",
      "and RT/DDRi/ICI cohorts; the monotherapy DDRi / ICI cohorts and",
      "cohorts 4-8 of DDRIO1835 were excluded. The published dataset (Data",
      "S1) has 212 mice and 2143 tumor-diameter records: control 57, RT 36,",
      "RT/PARPi 24, RT/ATMi 17, RT/ICI 42, RT/PARPi/ICI 24, RT/ATMi/ICI 12.",
      "Mice reaching a complete response were censored afterwards."
    )
  )

  ini({
    # --- Estimated structural parameters (Table 2, model 4) ---
    lrbase_core <- log(3.83); label("Initial quiescent-core depth Q0 (mm)") # Table 2 Q0 = 3.83 (RSE 2%)
    lp <- log(0.239); label("Baseline logistic rim growth rate lambda (1/day)") # Table 2 lambda = 0.239 (RSE 3%)
    lalpha_rt <- log(3.06); label("APC recruitment by RT alpha (cells/(Gy*mm*day) per mm rim)") # Table 2 alpha = 3.06 (RSE 2%)
    lkexh <- log(0.278); label("Natural CD8 T-cell exhaustion rate beta (1/day)") # Table 2 beta = 0.278 (RSE 23%)

    # Drug effects: each indicator adds gamma to the divisor of the exhaustion
    # rate, beta / (1 + gamma_PARPi + gamma_ATMi + gamma_ICI) (model 4).
    e_trt_olaparib_kexh <- 2.03; label("Olaparib (PARPi) inhibition of T-cell exhaustion gamma_PARPi (unitless)") # Table 2 gamma PARPi = 2.03 (RSE 46%)
    e_trt_azd0156_kexh <- 4.34; label("AZD0156 (ATMi) inhibition of T-cell exhaustion gamma_ATMi (unitless)") # Table 2 gamma ATMi = 4.34 (RSE 40%)
    e_trt_antipdl1_kexh <- 11.9; label("Anti-PD-L1 (ICI) inhibition of T-cell exhaustion gamma_ICI (unitless)") # Table 2 gamma ICI = 11.9 (RSE 32%)

    # --- Literature / assumed constants (Table 2, not estimated) ---
    rbase_rim <- fixed(2.4); label("Initial proliferating-rim depth P0 (mm)") # Table 2 P0 = 2.4 (ref 36); Code S1 A_0(1) = 2.4
    kcap_rim <- fixed(2.4); label("Rim carrying capacity K (mm)") # Table 2 K = 2.4 (ref 36)
    kpt <- fixed(0.001); label("CD8 killing rate in the rim KPT (1/(cell*day))") # Table 2 KPT = 0.001 (ref 21)
    kptq <- fixed(0.001); label("CD8-driven transfer of core into rim KPTQ (1/(cell*day))") # Table 2 KPTQ = 0.001 (ref 21)
    kr <- fixed(0.1); label("Affinity of core-to-rim transfer KR (mm)") # Table 2 KR = 0.1 (assumed, footnote b)
    kin_apc_cells <- fixed(0.001); label("Baseline APC influx K_phiA (cells/(mm*day))") # Table 2 K_phiA = 0.001 (assumed, footnote a)
    kout_apc_cells <- fixed(0.648); label("APC natural death rate K_Aphi (1/day)") # Table 2 K_Aphi = 0.648 (ref 37); Code S1 KDTL = 0.027*24
    kact_cd8 <- fixed(9.12); label("APC activation of CD8 T cells K_AT (1/day)") # Table 2 K_AT = 9.12 (ref 37); Code S1 KDLCD8 = 0.38*24

    # RT input of Equation 2, R(t) = 2 Gy * Heaviside(4.1 - t): the daily RT
    # dose and the time the Heaviside switches off. Experimental-design
    # constants, overridable to simulate other schedules.
    rt_dose <- fixed(2); label("Daily RT dose (Gy)") # Equation 2 (2 Gy); Table 1 footnote a
    # Code S1 $DES switches RT on with 'IF (TIME.LT.4)'. NONMEM's TIME there is
    # the time of the record being integrated towards, so on the Monday /
    # Wednesday / Friday measurement records (days 0, 2, 4, ...) RT is on only
    # for days 0-2: the window under which Table 2 was estimated, and the one
    # that reproduces the Figure 3 VPC and the Data S1 medians. Figure 4 was
    # simulated on a dense day 0-34 grid, where the same code gives a 4-day
    # window; set rt_end = 4 (or 4.1, Equation 2) to reproduce it.
    rt_end <- fixed(2); label("End of the RT Heaviside window (day)") # Code S1 'TIME.LT.4' on Data S1 record times 0, 2, 4

    # --- IIV: full 3x3 block on log Q0, log lambda, log beta (Table 2) ---
    # Variances 0.0429 / 0.064 / 1.25; correlations -0.27 (Q0,lambda),
    # -0.04 (Q0,beta), -0.11 (lambda,beta), each written as the covariance
    # r * sqrt(var_i * var_j).
    etalrbase_core + etalp + etalkexh ~ c(
      0.0429,
      -0.27 * sqrt(0.0429 * 0.064), 0.064,
      -0.04 * sqrt(0.0429 * 1.25), -0.11 * sqrt(0.064 * 1.25), 1.25
    )

    # --- Residual error ---
    # Table 2 'Additive RV' = 0.506 mm^2 is the NONMEM $SIGMA variance; SD = sqrt(0.506)
    addSd <- 0.71134; label("Additive residual error SD on tumor diameter (mm)") # Table 2 sigma = 0.506 mm^2 (RSE 4%)
  })

  model({
    # Individual parameters
    q0 <- exp(lrbase_core + etalrbase_core)
    lambda <- exp(lp + etalp)
    alpha_rt <- exp(lalpha_rt)
    kexh <- exp(lkexh + etalkexh)

    # RT input (Equation 2): rt_dose Gy for t < rt_end in RT-treated mice
    rt <- TRT_RADIOTHERAPY * rt_dose * (t < rt_end)

    # Exhaustion rate inhibited by the combination drugs (model 4)
    kexh_trt <- kexh /
      (1 + e_trt_olaparib_kexh * TRT_OLAPARIB + e_trt_azd0156_kexh * TRT_AZD0156 + e_trt_antipdl1_kexh * TRT_ANTIPDL1)

    # Equation 1: logistic proliferating rim P and quiescent core Q
    fq <- core / (kr + core)
    d/dt(rim) <- lambda * rim * (1 - rim / kcap_rim) + kptq * rim * cd8_cells * fq - kpt * rim * cd8_cells
    d/dt(core) <- lambda * rim^2 / kcap_rim - kptq * rim * cd8_cells * fq
    # Equation 2: APCs recruited by RT-induced cell death in the rim
    d/dt(apc_cells) <- kin_apc_cells + alpha_rt * rim * rt - kout_apc_cells * apc_cells
    # Model 4: CD8 T cells activated by APCs and exhausted at kexh_trt
    d/dt(cd8_cells) <- kact_cd8 * apc_cells - kexh_trt * cd8_cells

    # Initial conditions (Table 2): A0 = C0 = 0
    rim(0) <- rbase_rim
    core(0) <- q0
    apc_cells(0) <- 0
    cd8_cells(0) <- 0

    # Observation: tumor diameter = rim depth + core depth (Code S1 $ERROR)
    tumor_size <- rim + core
    tumor_size ~ add(addSd)
  })
}
